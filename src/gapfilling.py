import pandas as pd
from cobra.flux_analysis import gapfilling
import cobra
import os
from pathlib import Path
import subprocess
import tempfile
import gc
from Bio import SeqIO
import findrefnames
import pickle
from math import isfinite
from src.solver_config import ensure_solver
from src.gpr_transfer import replace_gene_token



def read_gapgenes(model, genes, solution, refname, name):
    findrefnames.predata(refname)
    output = Path('working') / name / f'outfile{name}.csv'
    df = pd.read_csv(output, sep='\t', names=range(12), comment='#')
    scores = pd.to_numeric(df[11], errors='raise')
    if not df.empty and (df.isna().any().any() or not scores.map(isfinite).all()):
        raise ValueError(f'Malformed gapfilling alignment output: {output}')
    df[11] = scores
    df = df[df[11] >= 200].sort_values(11, ascending=False, kind='stable')
    best_hits = df.drop_duplicates(0).set_index(0)[1].to_dict()
    nonfindgenes = []
    for gene in genes:
        accession = gene.split('|')[1] if '|' in gene else gene
        reference_gene = findrefnames.findentryname(accession)
        if not isinstance(reference_gene, str) or not reference_gene.strip():
            raise ValueError(f'Reference protein missing ordered locus mapping: {accession}')
        if gene not in best_hits:
            nonfindgenes.append(reference_gene)
            continue
        target = best_hits[gene]
        new_gene = target.split('|')[1] if '|' in target else target
        for reaction_id in solution:
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.gene_reaction_rule = replace_gene_token(
                reaction.gene_reaction_rule, reference_gene, new_gene,
            )
    cobra.manipulation.remove_genes(model, [gene for gene in model.genes if not gene.reactions])
    return model, nonfindgenes


def run_gapgene_alignment(name, refname):
    """Publish a new alignment only after the external commands succeed."""
    work = Path('working') / name
    database = work / f'{name}db'
    query = work / f'{refname}_gap.fasta'
    output = work / f'outfile{name}.csv'
    nucleotide = work / f'{name}_non_anno.fasta'
    with tempfile.TemporaryDirectory(dir=work, prefix='gapgene-') as temporary:
        pending = Path(temporary) / 'hits.tsv'
        if nucleotide.exists():
            commands = [
                ['makeblastdb', '-in', str(nucleotide), '-dbtype', 'nucl',
                 '-input_type', 'fasta', '-out', str(database)],
                ['tblastn', '-query', str(query), '-db', str(database), '-out', str(pending),
                 '-outfmt', '7', '-evalue', '1e-1', '-sorthits', '1'],
            ]
        else:
            commands = [
                ['diamond', 'makedb', '--in', str(work / f'{name}.fasta'), '--db', str(database)],
                ['diamond', 'blastp', '-q', str(query), '-d', str(database),
                 '--out', str(pending)],
            ]
        for command in commands:
            subprocess.run(command, check=True)
        if not pending.is_file():
            raise RuntimeError('Gapfilling alignment succeeded without producing an output file')
        pending.replace(output)

def configure_gapfiller(gapfiller, time_limit, soft_mem_limit, mip_gap):
    """Apply the same numerical/time/gap limits to Gurobi and COPT."""
    configuration = gapfiller.model.solver.configuration
    # optlang's "auto" retries every non-optimal status, even a time limit,
    # with a fresh native solve. Return resource-limited incumbents directly.
    configuration.presolve = True
    tolerances = configuration.tolerances
    tolerances.feasibility = 1e-6
    tolerances.integrality = 1e-6
    tolerances.optimality = 1e-6

    solver_name = gapfiller.model.solver.interface.__name__.lower()
    problem = gapfiller.model.solver.problem
    if 'gurobi' in solver_name:
        problem.setParam('TimeLimit', time_limit)
        problem.setParam('SoftMemLimit', soft_mem_limit)
        problem.setParam('MIPGap', mip_gap)
    elif 'copt' in solver_name:
        from coptpy import COPT

        problem.setParam(COPT.Param.TimeLimit, time_limit)
        # COPT MemLimit is in MB and is a hard limit, unlike SoftMemLimit.
        problem.setParam(COPT.Param.MemLimit, soft_mem_limit * 1024)
        problem.setParam(COPT.Param.RelGap, mip_gap)
    else:
        raise RuntimeError(f'Unsupported AlphaGEM gapfilling solver: {solver_name}')


def fill_with_validation_retry(gapfiller):
    """Retry only validation failures with stricter numerical tolerances."""
    try:
        return gapfiller.fill()[0]
    except RuntimeError as error:
        if ('Failed to validate gap filled model' not in str(error)
                or gapfiller.model.solver.status != 'optimal'):
            raise

        print('Gapfilling validation failed at 1e-6; retrying at 1e-9.')
        tolerances = gapfiller.model.solver.configuration.tolerances
        tolerances.feasibility = 1e-9
        tolerances.integrality = 1e-9
        tolerances.optimality = 1e-9
        gapfiller.integer_threshold = 1e-9
        return gapfiller.fill()[0]


def has_feasible_incumbent(model):
    """Require a native primal solution, including at a resource limit."""
    solver = model.solver
    if solver.status in {None, 'loaded', 'undefined', 'infeasible', 'nofeasible',
                         'unbounded', 'infeasible_or_unbounded', 'numeric',
                         'in_progress'}:
        return False
    native = solver.problem
    solver_name = solver.interface.__name__.lower()
    if 'gurobi' in solver_name:
        available = native.SolCount > 0
    elif 'copt' in solver_name:
        available = native.hasmipsol if solver.is_integer else native.haslpsol
    else:
        return False
    if not available:
        return False
    value = solver.objective.value
    return value is not None and isfinite(value)


def has_growth(model, solution, threshold=0):
    value = solution.objective_value
    return (value is not None and isfinite(value) and value > threshold
            and has_feasible_incumbent(model))


def require_growth(model, stage):
    solution = model.optimize()
    if not has_growth(model, solution):
        raise RuntimeError(
            f'{stage} has no finite positive growth from a feasible solution '
            f'(status={solution.status}, objective={solution.objective_value}).'
        )
    return solution.objective_value


def fill_with_validated_incumbent(gapfiller):
    """Accept a resource-limited MIP only after validating its incumbent."""
    try:
        return fill_with_validation_retry(gapfiller)
    except (cobra.exceptions.OptimizationError, RuntimeError) as error:
        if (isinstance(error, RuntimeError)
                and 'Failed to validate gap filled model' not in str(error)):
            raise
        status = gapfiller.model.solver.status
        if not has_feasible_incumbent(gapfiller.model):
            raise RuntimeError(
                f'Gapfilling failed without a feasible incumbent (status={status}).'
            ) from error
        primals = gapfiller.model.solver.primal_values.values()
        if any(value is None or not isfinite(value) for value in primals):
            raise RuntimeError('Gapfilling incumbent contains missing flux values.') from error
        reactions = []
        for indicator in gapfiller.indicators:
            value = indicator._get_primal()
            if value is None or not isfinite(value):
                raise RuntimeError('Gapfilling incumbent contains missing indicators.') from error
            if value > gapfiller.integer_threshold:
                reactions.append(gapfiller.model.reactions.get_by_id(indicator.rxn_id))
        try:
            valid = gapfiller.validate(reactions)
        except Exception as validation_error:
            raise RuntimeError('Could not validate the gapfilling incumbent.') from validation_error
        if not valid:
            raise RuntimeError(
                f'Gapfilling incumbent failed growth validation (status={status}).'
            ) from error
        print(f'Using validated gapfilling incumbent (status={status}).')
        return reactions


def gapfill(name, refname, refmodelname, grothmedium='min', esstinal=None,
            use_full=False, time_limit=3600, soft_mem_limit=8, mip_gap=0.01):
    ensure_solver()
    if esstinal is None:
        esstinal=[]
    refmodel=cobra.io.read_sbml_model(f'models/{refmodelname}')
    # reactionpool=cobra.io.load_json_model('./models/reaction_pool_final.json')
    findrefnames.predata(refname)
    model = cobra.io.load_yaml_model(f'working/{name}/{name}-GEM_withgaps.yml')
    if refname == 'human':
        model.objective=model.reactions.MAR00021
        model.reactions.MAR00021.bounds=[0,1000]
    try:
        if has_growth(model, model.optimize(), threshold=0.001):
            cobra.io.write_sbml_model(model, f'./working/{name}/{name}-GEM.xml')
            cobra.io.save_json_model(model, f'./working/{name}/{name}-GEM.json')
            return 0
    except cobra.exceptions.OptimizationError:
        pass
    # gap = gapfilling.GapFiller(model, universal=refmodel, integer_threshold=1e-12, demand_reactions=False,lower_bound=0.05)
    # for j in gap.indicators:
    #     i=gap.model.variables.get(key=j.rxn_id)
    #     if i._get_primal() >0:
    #         solution.append(i.name)
    # refmodel_reactions=[rrr.id for rrr in refmodel.reactions]
    # for reac in model.reactions:
    #     if reac.id not in refmodel_reactions:
    #       refmodel.add_reactions([reac])
    # fba=cobra.flux_analysis.pfba(refmodel)
    # model2 = cobra.io.load_yaml_model(f'working/{name}/{name}-GEM_withgaps.yml')
    model2=model.copy()
    fastafile=open(f'./working/{name}/{refname}_gap.fasta','w')#要保存的gapgene的fasta文件
    targenes=[]
    allgenes=[x for x in SeqIO.parse(open(f'data_available/{refname}.fasta'),'fasta')]#打开fasta文件
    allgenesid=[x.id.split('|')[1] for x in allgenes ]
    reacgaps=[]
    # for solution1 in solution:
    #     if fba.fluxes.get(solution1)>=1e-10 or fba.fluxes.get(solution1)<=-1e-10:
    #         reacgaps.append(gap.model.reactions.get_by_id(solution1))
    objectbound=require_growth(refmodel, 'Reference model')
    gapfiller = gapfilling.GapFiller(model2, refmodel, integer_threshold=1e-6, demand_reactions=False,lower_bound=min(objectbound*0.8,0.09))
    configure_gapfiller(gapfiller, time_limit, soft_mem_limit, mip_gap)
    reacgaps=fill_with_validated_incumbent(gapfiller)
    for r in esstinal:
        try:
            reacgaps.append(refmodel.reactions.get_by_id(r))
        except:
            print('Can not find essential reaction {r}',r)
            continue
    for reac in reacgaps:
        for genes in reac.genes:
            try:
              targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
            except:
                continue
    medium = refmodel.medium
    if refname!='human' and refname!='strco' and refname!='synechocystis' and use_full==True:
        # GapFiller retains a large copied MIP model. Release the first solve
        # before building the full-medium MIP to avoid a doubled memory peak.
        del gapfiller
        gc.collect()
        with open(f'data_available/{refname}_full_medium.pkl', 'rb') as file:
            fullmedium = pickle.load(file)
        medium=refmodel.medium.copy()
        refmodel.medium = fullmedium
        objectbound = require_growth(refmodel, 'Full-medium reference model')
        fba2=cobra.flux_analysis.pfba(refmodel)
        reacgaps2 = []
        model2.medium=fullmedium
        # for solution1 in solution:
        #     if fba2.fluxes.get(solution1) >= 1e-10 or fba2.fluxes.get(solution1) <= -1e-10:
        #         reacgaps2.append(gap.model.reactions.get_by_id(solution1))
        gapfiller = gapfilling.GapFiller(model2, refmodel, integer_threshold=1e-6, demand_reactions=False,lower_bound=min(objectbound*0.5,0.05))
        configure_gapfiller(gapfiller, time_limit, soft_mem_limit, mip_gap)
        reacgaps2 = fill_with_validated_incumbent(gapfiller)
        if refname=='yeast':
           reacgaps2.append(refmodel.reactions.get_by_id('rhea_21624_c'))
        for reac in reacgaps2:
            for genes in reac.genes:
                try:
                    targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
                except:
                    continue
        #pd.DataFrame(reacgaps).to_excel(f'working/{name}/reacgapsMM.xlsx')
        #pd.DataFrame(reacgaps2).to_excel(f'working/{name}/reacgaps2FULL.xlsx')
        reacgaps=reacgaps+reacgaps2
    reacgaps=[refmodel.reactions.get_by_id(r) for r in list(set([r.id for r in reacgaps]))]
    model2.add_reactions(reacgaps)
    solution = model2.optimize()
    print(solution)
    if not has_growth(model2, solution, threshold=1e-5):
        gapnow=[]
        for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            with model2:
                model2.add_reactions([i])
                if has_growth(model2, model2.optimize(), threshold=1e-4):
                    gapnow.append(i)
        # print(gapnow)
        for reac in gapnow:
            for genes in reac.genes:
                try:
                    targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
                except:
                    continue
        model2.add_reactions(gapnow)
        print(f'Adding {len(gapnow)} genes in gapfilling')
    nonfindgenes=[]
    if len(targenes)!=0:
       SeqIO.write(targenes,fastafile,'fasta')
       fastafile.close()
       run_gapgene_alignment(name, refname)
       model2,nonfindgenes=read_gapgenes(model=model2,genes=[x.id for x in targenes],solution=[reac.id for reac in reacgaps],refname=refname,name=name)
    model2.optimize()
    # print([reac.id for reac in reacgaps])
       # for met in model.metabolites:
       #     if met.compartment == '':
       #         model.metabolites.get_by_id(met.id).compartment = 'c'
       # for re in model.reactions:
       #     for cp in re.compartments:
       #         if cp == '':
       #             model.reactions.get_by_id(re.id).compartments.remove('')
       #             model.reactions.get_by_id(re.id).compartments.add('c')
    if grothmedium=='min':
        model2.medium=medium
    #changed
    # with model2:
    #     nonfindreactions=[r.id for r in cobra.manipulation.delete.knock_out_model_genes(model2, nonfindgenes)]
    #     pd.DataFrame(nonfindreactions).to_excel(f'working/{name}/reacgaps_non_gene.xlsx')
    filter1=min(require_growth(model2, 'Gapfilled model before gene filtering'),0.1)
    num_nonreactions=[]
    for gene in nonfindgenes:
        try:
            geneid=model2.genes.get_by_id(gene)
        except:
            continue
        with model2:
            cobra.manipulation.delete.knock_out_model_genes(model2,[geneid])
            try:
              if not has_growth(model2, model2.optimize(), threshold=0.5*filter1):
                num_nonreactions=num_nonreactions+[r.id for r in geneid.reactions]
                continue
            except:
                num_nonreactions=num_nonreactions+[r.id for r in geneid.reactions]
                continue
        # print('delete gene ',gene)
        cobra.manipulation.delete.remove_genes(model2,[geneid])
    print(model2.optimize().objective_value)
    model2.compartments=refmodel.compartments
    for r in esstinal:
        try:
            model2.add_reactions([refmodel.reactions.get_by_id(r)])
        except:
            print('Can not find essential reaction {r}',r)
            continue
    require_growth(model2, 'Final gapfilled model')
    # if refname=='ecoli':
    #     cobra.io.save_json_model(model2,f'./working/{name}/{name}-GEM.json')
    #     model2 = cobra.io.load_json_model(f'./working/{name}/{name}-GEM.json')
    #     model2.id = f'{name}'
    #     model2.name = f'{name}-GEM'
    #     for i in model2.reactions:
    #         keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
    #         # 删除这些键
    #         for k in keys_to_remove:
    #             del i.annotation[k]
    #     cobra.io.save_json_model(model2,f'./working/{name}/{name}-GEM.json')
    #     cobra.io.write_sbml_model(model2, f'./working/{name}/{name}-GEM.xml')
    # else:
    for i in model2.reactions:
        keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
        # 删除这些键
        for k in keys_to_remove:
                del i.annotation[k]
    cobra.io.save_json_model(model2, f'./working/{name}/{name}-GEM.json')
    cobra.io.write_sbml_model(model2, f'./working/{name}/{name}-GEM.xml')
    model2 = cobra.io.read_sbml_model(f'./working/{name}/{name}-GEM.xml')
    require_growth(model2, 'Reloaded final model')
    model2.id = f'{name}'
    model2.name = f'{name}-GEM'
    for gene_anno in model2.genes:
        gene_anno.annotation['uniprot']=gene_anno.id
    cobra.io.write_sbml_model(model2, f'./working/{name}/{name}-GEM.xml')
    cobra.io.save_json_model(model2, f'./working/{name}/{name}-GEM.json')

    print(f'Number of all gapfilling reactions: {len(reacgaps)}')
    print(f'Number of non-bio genes: {len(set([r.id for r in reacgaps])&set(num_nonreactions))}')

