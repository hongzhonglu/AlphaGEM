import cobra
import pandas as pd
def datapre(refmodel,refname):
    global refer
    global transgenes
    refer = pd.read_excel(f'data_available/{refname}.xlsx')
    if refname == 'yeast' or refname == 'strco' or refname == 'human':
        ymodel = cobra.io.read_sbml_model(f'models/{refmodel}')
    if refname == 'ecoli':
        ymodel = cobra.io.load_json_model('models/iML1515.json')
    transgenes = []
    for reaction in ymodel.reactions:
        reac = []
        pro=[]
        for reactant in reaction.reactants:
            reac.append(reactant)
        for product in reaction.products:
            pro.append(product)
        #for meta in reaction.metabolites:
            #reac.append(meta.name)
        #if 2 * len(set(reac)) == len(reaction.metabolites):
            #for gene in reaction.genes: transgenes.append(gene.id)
            #continue
        try:
            if len(reac)==len(pro)==1 and reac[0].name==pro[0].name and reac[0].compartments!=pro[0].compartments:
                for gene in reaction.genes: transgenes.append(gene.id)
                continue
        except:
            0
        if 'transport' in reaction.name or 'permease' in reaction.name or 'channel' in reaction.name or 'pore' in reaction.name or 'uniport' in reaction.name or 'symport' in reaction.name or 'antiport' in reaction.name or 'gradient' in reaction.name or 'ABC' in reaction.name or 'ATPase' in reaction.name or 'PTS' in reaction.name:
            for gene in reaction.genes: transgenes.append(gene.id)
    transgenes= list(set(transgenes))
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

