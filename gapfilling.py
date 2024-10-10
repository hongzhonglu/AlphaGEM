import pandas as pd
from cobra.flux_analysis import gapfilling
import cobra
import os
from Bio import SeqIO
import findrefnames
import model_build1
def read_gapgenes(model,genes,solution,refname,name):
    findrefnames.predata(refname)
    df=pd.read_csv(f'ziyuan/outfile{name}.csv',sep='\t',names=range(12))
    for gene in genes:
        try:
            newgenename=df.iat[list(df[0]).index(gene),1].split('|')[1]
            for solution1 in solution:
                try:
                   model.reactions.get_by_id(solution1).gene_reaction_rule=model.reactions.get_by_id(solution1).gene_reaction_rule.replace(findrefnames.findentryname(gene.split('|')[1]),newgenename)
                   print(model.reactions.get_by_id(solution1).gene_reaction_rule)
                except:
                    continue
        except:
            continue
    return model

def gapfill(name,refname):
    if refname=='yeast':
        refmodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
    if refname=='ecoli':
        refmodel=cobra.io.load_json_model('models/iML1515.json')
    if refname=='strco':
        refmodel=cobra.io.read_sbml_model('models/Sco-GEM.xml')
    if refname=='human':
        refmodel = cobra.io.read_sbml_model('models/Human-GEM.xml')
    solution=[]
    findrefnames.predata(refname)
    model = cobra.io.load_yaml_model(f'models/tarmodel_{name}.yml')
    if refname == 'human':
        model.objective=model.reactions.MAR00021
    gap = gapfilling.GapFiller(model, universal=refmodel, integer_threshold=1e-12, demand_reactions=False,lower_bound=0.05)
    for j in gap.indicators:
        i=gap.model.variables.get(key=j.rxn_id)
        if i._get_primal() >=0:
            solution.append(i.name)
    for solution1 in solution:
        model.add_reactions([gap.model.reactions.get_by_id(solution1)])
    fba = model.optimize()
    model = cobra.io.load_yaml_model(f'models/tarmodel_{name}.yml')
    fastafile=open(f'ziyuan/{refname}_gap.fasta','w')#要保存的gapgene的fasta文件
    targenes=[]
    allgenes=[x for x in SeqIO.parse(open(f'ziyuan/{refname}.fasta'),'fasta')]#打开fasta文件
    allgenesid=[x.id.split('|')[1] for x in allgenes ]
    reacgaps=[]
    for solution1 in solution:
        if fba.fluxes.get(solution1)>=1e-10 or fba.fluxes.get(solution1)<=-1e-10:
            reacgaps.append(gap.model.reactions.get_by_id(solution1))
    for reac in reacgaps:
        for genes in reac.genes:
            try:
              targenes.append(allgenes[allgenesid.index(findrefnames.findmodelname(genes.id))])
            except:
                continue
    model.add_reactions(reacgaps)
    model.optimize()
    if len(targenes)!=0:
       SeqIO.write(targenes,fastafile,'fasta')
       fastafile.close()
       os.system(f'diamond makedb --in ziyuan/{name}.fasta --db ziyuan/{name}db')
       os.system(f'diamond blastp -q ziyuan/{refname}_gap.fasta -d ziyuan/{name}db --out ziyuan/outfile{name}.csv')
    print([x.id for x in targenes])
    model=read_gapgenes(model=model,genes=[x.id for x in targenes],solution=[reac.id for reac in reacgaps],refname=refname,name=name)
    model.optimize()
    print([reac.id for reac in reacgaps])
    for met in model.metabolites:
        if met.compartment == '':
            model.metabolites.get_by_id(met.id).compartment = 'c'
    for re in model.reactions:
        for cp in re.compartments:
            if cp == '':
                model.reactions.get_by_id(re.id).compartments.remove('')
                model.reactions.get_by_id(re.id).compartments.add('c')
    model.optimize()
    cobra.io.write_sbml_model(model, f'models/tarmodel__{name}text.xml')




