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


def model_reaction(name,refmodel):
    refmodel=cobra.io.read_sbml_model(f'models/{refmodel}')
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
    metdict={metss.name.lower():metss for metss in model.metabolites if metss.compartment=='c'}
    for metss in model.metabolites:
        if metss.compartment=='c':
            metdict[metss.name.lower().replace('(','').replace(')','')]=metss
    for index, row in tqdm(gpr_e.iterrows()):
        reactions=[]
        reactions.append(str(row['rhea']))
        for rea in reactions:
            try:
                model.reactions.get_by_id('rhea_'+str(rea)+'_c')#已经添加过了
                continue
            except:
                pass
            already =0
            try:
               reaction_add=reactionpools.reactions.get_by_id('rhea_'+str(rea)+'_c')
            except:
                print('rhea_'+str(rea)+'_c'+' does not exist in reaction pool')
                continue
            if (reaction_add.reactants==[] or reaction_add.products==[]) and reaction_add.id not in refids:
                print('reaction error')
                continue
            if reaction_add.id in refids:
                reaction_add=refmodel.reactions.get_by_id(reaction_add.id)
                already=1
            reaction_add.gene_reaction_rule=row['reaction']
            if reaction_add.id=='rhea_24628_c' or reaction_add.id=='rhea_24632_c':
                continue
            reaction_add.annotation['from']='NonHomo'
            #TODO: DEBUG REACTION ERROR OCCUR
            a=list(reaction_add.metabolites.keys())
            for met in a:
                 if already==1:
                     continue
                 try:
                     model.metabolites.get_by_id(met)
                     continue
                 except:
                     pass
                 metsearch=met.name.lower() if met.name.lower() in list(metdict.keys()) else met.name.lower().replace('(','').replace(')','')
                 try:
                     if metdict[metsearch].elements==met.elements and metdict[
                         metsearch].charge==met.charge and metdict[metsearch].compartment==met.compartment:
                         stoich=reaction_add.metabolites[met]
                         reaction_add.add_metabolites({metdict[metsearch]:stoich})
                         reaction_add.add_metabolites({met:-1*stoich})
                         continue
                     # try:
                     #   if (metdict[metsearch].elements['H'] - met.elements['H'] == metdict[
                     #     metsearch].charge - met.charge and metdict[metsearch].elements.pop(
                     #     'H') == met.elements.pop('H')) and metdict[metsearch].compartment==met.compartment:
                     #     stoich = reaction_add.metabolites[met]
                     #     reaction_add.add_metabolites({metdict[metsearch]: reaction_add.metabolites[stoich],
                     #                                   met: -1 * reaction_add.metabolites[stoich]})
                     # except:
                     #   pass
                 except KeyError as e:
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