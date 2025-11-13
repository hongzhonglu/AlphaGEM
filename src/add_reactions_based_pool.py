import cobra
import pandas as pd
import pickle
from tqdm import tqdm
def check_reactions(reac,model):
    result={}
    for met in reac.metabolites.keys():
        metvalue=reac.metabolites[met]
        try:
            metabolite=model.metabolites.get_by_id(met).elements
        except:
            metabolite=met.elements
        result={k: result.get(k, 0) + metvalue*metabolite.get(k, 0) for k in set(result) | set(metabolite)}
    for i in result:
        if result[i]!=0:
            return 0
    return 1


# def reactionpoolcomp(compartment,transportout,reactionpool):
#     for reaction in reactionpool.reactions:
#         if len(reaction.compartments) == 1:
#             reaction.id=reaction.id+'_'+compartment
#         else:
#             reaction.id=reaction.id+'_'+compartment
#     for metbolites in reactionpool.metabolites:
#         if metbolites.compartment == 'out':
#             metbolites.id=metbolites.id.split('_')[0]+'_'+transportout
#             metbolites.compartment=transportout
#         else:
#             metbolites.id=metbolites.id.split('_')[0]+'_'+compartment
#             metbolites.compartment=compartment
#     return reactionpool



def model_reaction(name,refname):
    if refname=='yeast':
        refmodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
    if refname=='ecoli':
        refmodel=cobra.io.load_json_model('models/iML1515.json')
    if refname=='strco':
        refmodel=cobra.io.read_sbml_model('models/Sco-GEM.xml')
    if refname=='human':
        refmodel = cobra.io.read_sbml_model('models/Human-GEM.xml')
    refids=[r.id for r in refmodel.reactions]
    model = cobra.io.load_yaml_model(f'./working/{name}/model2{name}.yml')
    fba_begin=model.optimize().objective_value
    gpr_e = pd.read_excel(f'./working/{name}/{name}gpr_score.xlsx')
    pklopen=open('data_available/general2specific.pkl', 'rb')
    general2spe=pickle.load(pklopen)
    reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
    pklopen.close()
    Reactions=[]
    reactionpools=reactionpool.copy()
    # reactionpools=reactionpoolcomp('c','e',reactionpools)
    metdict={metss.name.lower():metss for metss in model.metabolites if metss.compartment=='c'}
    for metss in model.metabolites:
        if metss.compartment=='c':
            metdict[metss.name.lower().replace('(','').replace(')','')]=metss
    for index, row in tqdm(gpr_e.iterrows()):
        reactions=[]
        reactions.append(str(row['rhea']))
        while True:
            found = False
            for item in reactions[:]:
                if item in general2spe:
                    found = True
                    reactions.remove(item)
                    reactions.extend(general2spe[item])
            if not found:
                break
        for rea in reactions:
            reaction_add=reactionpools.reactions.get_by_id('rhea_'+str(rea)+'_c')
            if (reaction_add.reactants==[] or reaction_add.products==[]) and reaction_add.id not in refids:
                print('reaction error')
                continue
            if reaction_add.id in refids:
                reaction_add=refmodel.reactions.get_by_id(reaction_add.id)
            reaction_add.gene_reaction_rule=row['reaction']
            if reaction_add.id=='rhea_24628_c' or reaction_add.id=='rhea_24632_c':
                continue
            reaction_add.annotation['from']='NonHomo'
            if reaction_add.id in refids:
                reaction_add=refmodel.reactions.get_by_id(reaction_add.id)
            for met in reaction_add.metabolites.keys():
                 metsearch=met.name.lower() if met.name.lower() in list(metdict.keys()) else met.name.lower().replace('(','').replace(')','')
                 try:
                     if metdict[metsearch].elements==met.elements and metdict[
                         metsearch].charge==met.charge and metdict[metsearch].compartment==met.compartment:
                         reaction_add.add_metabolites({metdict[metsearch]:reaction_add.metabolites[met],
                                                       met:-1*reaction_add.metabolites[met]})
                         print(f'{metsearch} is added to {reaction_add.id}')
                         print(reaction_add.metabolites)
                         continue
                     try:
                       if (metdict[metsearch].elements['H'] - met.elements['H'] == metdict[
                         metsearch].charge - met.charge and metdict[metsearch].elements.pop(
                         'H') == met.elements.pop('H')) and metdict[metsearch].compartment==met.compartment:
                         reaction_add.add_metabolites({metdict[metsearch]: reaction_add.metabolites[met],
                                                       met: -1 * reaction_add.metabolites[met]})
                         print(f'{metsearch} is added to {reaction_add.id}')
                         print(reaction_add.metabolites)
                     except:
                       pass
                 except KeyError as e:
                     print(e)
                     pass
            try:
                if check_reactions(reaction_add, model)!=0:
                    model.add_reactions([reaction_add])
            except:
                pass
    model.id=f'{name}-GEM'
    for i in model.reactions:
        keys_to_remove = [k for k, v in i.annotation.items() if v in [[], '']]
        # 删除这些键
        for k in keys_to_remove:
            del i.annotation[k]
    cobra.io.save_yaml_model(model, f'./working/{name}/{name}-GEM_withgaps.yml')