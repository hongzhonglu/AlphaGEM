import cobra
import pandas as pd
from config import refmodel
from config import refname
refer=pd.read_excel(f'ziyuan/{refname}.xlsx')
if refname=='yeast' or refname== 'strco' or refname=='human':
    ymodel = cobra.io.read_sbml_model(f'models/{refmodel}')
if refname=='ecoli':
    ymodel=cobra.io.load_json_model('models/iML1515.json')
transgenes=[]
for reaction in ymodel.reactions:
    reac=[]
    for meta in reaction.metabolites:
         reac.append(meta.name)
    if 2*len(set(reac))==len(reaction.metabolites):
        for gene in reaction.genes: transgenes.append(gene.id)
        continue
    if 'transport' in reaction.name or 'aminotransferase' in reaction.name :
        for gene in reaction.genes:transgenes.append(gene.id)
def transpoter(genes):
    try:
        name = refer.iat[list(refer['Entry']).index(genes), 2]
        try:
            transgenes.index(name)
            return 1
        except:
            return 0
    except:
        print('error')

