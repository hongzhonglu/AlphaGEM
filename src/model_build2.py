import cobra
import pandas as pd
from tqdm import tqdm
def str_chanege_for_gpr(old_gpr,oldgene,newgenes):
    new_gpr=''
    oldgprcell=old_gpr.split(' or ')
    gprcell=[]
    for a in oldgprcell:
        if oldgene in a:
            for genes in newgenes:
                if 'and' in a:
                    gprcell.append('('+a.replace(oldgene,genes).replace('(','').replace(')','')+')')
                else:
                    gprcell.append(a.replace(oldgene, genes).replace('(', '').replace(')', ''))
        else:
            if 'and' in a:
                gprcell.append('('+a.replace('(','').replace(')','')+')')
            else:
                gprcell.append(a.replace('(', '').replace(')', ''))
    new_gpr=' or '.join(list(set(gprcell)))
    return new_gpr


def get_complex(leftreactions,homodict):
    complexdict={}
    for reaction in tqdm(leftreactions,'preparing other reactions'):
        for gene in list(reaction.genes):
            try:
                reaction.gene_reaction_rule=str_chanege_for_gpr(reaction.gene_reaction_rule,gene.id,homodict[gene.id])
            except:
                reaction.gene_reaction_rule=str_chanege_for_gpr(reaction.gene_reaction_rule,gene.id,['None'])
        complexdict[reaction.id]=[gpr.replace('(','').replace(')','') for gpr in reaction.gene_reaction_rule.split(' or ')]
    return complexdict

def modelbuild(refmodel,name):
    tarname = name
    tarmod = cobra.io.load_yaml_model(f'working/{name}/model1{name}.yml')
    print('imported target model')
    if refmodel=='yeast-GEM.xml' or refmodel=='Sco-GEM.xml' or refmodel=='Human-GEM.xml':
      ymod = cobra.io.read_sbml_model(f'models/{refmodel}')
    if refmodel=='iML1515.json':
        ymod=cobra.io.load_json_model(f'models/{refmodel}')
    tarmodel = cobra.Model('Target_model')
    homos = pd.read_excel(f'working/{name}/matrix_homolog{name}.xlsx')
    addedreactions=[reaction.id for reaction in tarmod.reactions]
    leftreactions=[]
    complexdict={}
    for reaction in tqdm(ymod.reactions,'finding left reactions'):
        if reaction.id in addedreactions:
            continue
        if 'and' in reaction.gene_reaction_rule:
            leftreactions.append(reaction)
    homos.columns = ['number', 'refmodelgene', 'tarmodelgene']
    homodict={}
    for groups in homos.groupby('refmodelgene'):
        homodict[groups[0]]=list(groups[1]['tarmodelgene'])
    #add othergenes
    for i in tqdm(range(len(tarmod.reactions))):
        tarmodel.add_reactions([tarmod.reactions[i]])

    complexdict=get_complex(leftreactions,homodict)
    complextruedict={}
    for key,gpr in tqdm(complexdict.items(),'adding complex'):
        complextruedict[key]=[]
        for gprs in gpr:
            if gprs.split(' and ').count('None')<=0.2*len(gprs.split(' and ')):
                complextruedict[key].append('('+' and '.join([x for x in gprs.split(' and ') if x != 'None'])+')')
        if complextruedict[key]==[]:
            complextruedict.pop(key)
    toaddreactions=[]
    for key in tqdm(complextruedict.keys(),'changing gpr'):
        if len(complextruedict[key])>1:
            ymod.reactions.get_by_id(key).gene_reaction_rule=' or '.join(complextruedict[key])
        else:
            ymod.reactions.get_by_id(key).gene_reaction_rule=complextruedict[key][0][1:-1]
        toaddreactions.append(ymod.reactions.get_by_id(key))
    print(f'added complexes {len(toaddreactions)}')
    tarmodel.add_reactions(toaddreactions)
    tarmodel.objective=ymod.objective
    fba = tarmodel.optimize()
    cobra.io.save_yaml_model(tarmodel, f'working/{tarname}/model2{tarname}.yml')
