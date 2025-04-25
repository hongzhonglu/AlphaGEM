import pandas as pd
from cobra.flux_analysis import gapfilling
import cobra
import os
from Bio import SeqIO
import findrefnames
import pickle



def read_gapgenes(model,genes,solution,refname,name):
    findrefnames.predata(refname)
    df=pd.read_csv(f'working/{name}/outfile{name}.csv',sep='\t',names=range(12))
    nonfindgenes=[]
    for gene in genes:
        try:
            newgenename=df.iat[list(df[0]).index(gene),1].split('|')[1]
            for solution1 in solution:
                try:
                   model.reactions.get_by_id(solution1).gene_reaction_rule=model.reactions.get_by_id(solution1).gene_reaction_rule.replace(findrefnames.findentryname(gene.split('|')[1]),newgenename)
                except:
                    continue
        except:
            nonfindgenes.append(findrefnames.findentryname(gene.split('|')[1]))
    return model,nonfindgenes

def gapfill(name,refname,grothmedium='min'):
    if refname=='yeast':
        refmodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
    if refname=='ecoli':
        refmodel=cobra.io.load_json_model('models/iML1515.json')
    if refname=='strco':
        refmodel=cobra.io.read_sbml_model('models/Sco-GEM.xml')
    if refname=='human':
        refmodel = cobra.io.read_sbml_model('models/Human-GEM.xml')
    findrefnames.predata(refname)
    model = cobra.io.load_yaml_model(f'working/{name}/{name}-GEM_withgaps.yml')
    if refname == 'human':
        model.objective=model.reactions.MAR00021
        model.reactions.MAR00021.bounds=[0,1000]
    try:
      if model.optimize().objective_value>0.001:
        cobra.io.write_sbml_model(model, f'./working/{name}/{name}-GEM.xml')
        return 0
    except:
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
    gapfiller = gapfilling.GapFiller(model, refmodel, integer_threshold=1e-09, demand_reactions=False)
    gapfiller.model.solver.configuration.tolerances.feasibility = 1e-09
    gapfiller.model.solver.configuration.tolerances.integrality = 1e-09
    gapfiller.model.solver.configuration.tolerances.optimality = 1e-09
    try:
        reacgaps= gapfiller.fill()[0]
    except:
        reacgaps=[]
        fba=gapfiller.model.optimize()
        for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            if fba.fluxes[i.id] != 0:
                try:
                    reacgaps.append(refmodel.reactions.get_by_id(i.id))
                except:
                    continue
    for reac in reacgaps:
        for genes in reac.genes:
            try:
              targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
            except:
                continue
    medium = refmodel.medium
    if refname!='human' and refname!='strco':
        with open(f'data_available/{refname}_full_medium.pkl', 'rb') as file:
            fullmedium = pickle.load(file)
        medium=refmodel.medium
        refmodel.medium = fullmedium
        fba2=cobra.flux_analysis.pfba(refmodel)
        reacgaps2 = []
        model2.medium=fullmedium
        # for solution1 in solution:
        #     if fba2.fluxes.get(solution1) >= 1e-10 or fba2.fluxes.get(solution1) <= -1e-10:
        #         reacgaps2.append(gap.model.reactions.get_by_id(solution1))
        gapfiller = gapfilling.GapFiller(model, refmodel, integer_threshold=1e-09, demand_reactions=False)
        gapfiller.model.solver.configuration.tolerances.feasibility = 1e-09
        gapfiller.model.solver.configuration.tolerances.integrality = 1e-09
        gapfiller.model.solver.configuration.tolerances.optimality = 1e-09
        try:
            reacgaps2 = gapfiller.fill()[0]
        except:
          reacgaps2 = []
          fba = gapfiller.model.optimize()
          for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            if fba.fluxes[i.id] != 0:
                try:
                    reacgaps2.append(refmodel.reactions.get_by_id(i.id))
                except:
                    continue
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
    print(model2.optimize())
    if model2.optimize().objective_value<=1e-5:
        gapnow=[]
        for j in gapfiller.indicators:
            i = gapfiller.model.reactions.get_by_id(j.rxn_id)
            with model2:
                model2.add_reactions([i])
                if model2.optimize().objective_value>1e-4:
                    gapnow.append(i)
        print(gapnow)
        for reac in gapnow:
            for genes in reac.genes:
                try:
                    targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
                except:
                    continue
        model2.add_reactions(gapnow)
    nonfindgenes=[]
    if len(targenes)!=0:
       SeqIO.write(targenes,fastafile,'fasta')
       fastafile.close()
       if os.path.exists(f'working/{name}/{name}_non_anno.fasta'):
           os.system(f'makeblastdb -in working/{name}/{name}_non_anno.fasta -dbtype nucl -input_type fasta -out working/{name}/{name}db')
           print(f'tblastn -query working/{name}/{refname}_gap.fasta -db working/{name}/{name}db -out working/{name}/outfile{name}.csv -outfmt 7 -evalue 1e-3')
           os.system(f'tblastn -query working/{name}/{refname}_gap.fasta -db working/{name}/{name}db -out working/{name}/outfile{name}.csv -outfmt 7 -evalue 1e-3')
       else:
           os.system(f'diamond makedb --in working/{name}/{name}.fasta --db working/{name}/{name}db')
           os.system(f'diamond blastp -q working/{name}/{refname}_gap.fasta -d working/{name}/{name}db --out working/{name}/outfile{name}.csv')
       model2,nonfindgenes=read_gapgenes(model=model2,genes=[x.id for x in targenes],solution=[reac.id for reac in reacgaps],refname=refname,name=name)
    model2.optimize()
    print([reac.id for reac in reacgaps])
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
    filter1=model2.optimize().objective_value#changed
    # with model2:
    #     nonfindreactions=[r.id for r in cobra.manipulation.delete.knock_out_model_genes(model2, nonfindgenes)]
    #     pd.DataFrame(nonfindreactions).to_excel(f'working/{name}/reacgaps_non_gene.xlsx')
    for gene in nonfindgenes:
        try:
            geneid=model2.genes.get_by_id(gene)
        except:
            continue
        with model2:
            cobra.manipulation.delete.knock_out_model_genes(model2,[geneid])
            try:
              if model2.optimize().objective_value<=0.2*filter1:#changed
                continue
            except:
                continue
        print('delete gene ',gene)
        cobra.manipulation.delete.remove_genes(model2,[geneid])
    print(model2.optimize().objective_value)
    if refname=='ecoli':
        cobra.io.save_json_model(model2,f'./working/{name}/{name}-GEM.json')
    cobra.io.write_sbml_model(model2, f'./working/{name}/{name}-GEM.xml')




