import pandas as pd
from bioservices import KEGG
import cobra
from cobra import io
import re
import numpy as np
import time
import random
cpd=pd.read_excel('ziyuan/cpd.xlsx')
def check_reactions(reac,model):
    result={}
    for met in reac:
        metvalue=reac[met]
        metabolite=model.metabolites.get_by_id(met).elements
        result={k: result.get(k, 0) + metvalue*metabolite.get(k, 0) for k in set(result) | set(metabolite)}
    print(result)
    for i in result:
        if result[i]!=0:
            return 0
    return 1

def find_tarname(name,strr):
    global cpd
    if np.sum(cpd['ID'].str.find(strr)!=-1)!=0:
        if name=='CHARGE':
            return eval(cpd[name][list(cpd['ID'].str.find(strr)).index(0)])
        return cpd[name][list(cpd['ID'].str.find(strr)>0).index(0)]
    else:
        if strr=='C00001':
            return 'H2O'
        if name=='CHARGE':
            return 0
        reponse = K.get('{}'.format(strr))
        met = reponse.split('\n')
        for i in met:
            if name in i:
                return i[12:i.find(';')]
    return 'none'
def find_tarname2(name,strr):
    global cpd
    if strr == 'C00001':
        return 'H2O'
    if np.sum(cpd['ID'].str.find(strr)!=-1)!=0:
        return cpd[name][list(cpd['ID'].str.find(strr)).index(0)]
    else:
        reponse = K.get('{}'.format(strr))
        met = reponse.split('\n')
        for i in met:
           if name in i:
               if i.find('(')==-1:
                   return [i[12:],'no']
               return [i[12:i.find('(')],'yes']
    return 'none'
def modelreaction(name):
    global K
    K = KEGG()
    global cpd
    model = io.load_yaml_model(f'models/tarmodel{name}.yml')
    reactions = pd.read_excel(f'juzhen/reactions_kegg{name}.xlsx')
    gpr = pd.read_excel(f'juzhen/gpr{name}.xlsx')
    model_reaction = []
    for i in model.metabolites:
        model_reaction.append(i.name)
    for index, row in reactions.iterrows():
        pan = 0
        time.sleep(random.random())
        Reaction = cobra.Reaction(row['REACTION'])
        Reaction.name = row['NAME']
        Reaction.id = row['REACTION']
        if type(row['EQUATION']) == float:
            continue
        Reaction.lower_bound = -1000
        Reaction.upper_bound = 1000
        reac_metabolic = {}
        metbs = row['EQUATION'].split('<=>')
        metb_left = metbs[0].split('+')
        metb_right = metbs[1].split('+')
        for metb in metb_left:
            model_reaction = []
            for i in model.metabolites:
                model_reaction.append(i.name.lower())
            time.sleep(random.random() * 5)
            metb = metb.replace('(n-2)', '2').replace('(m)', '').replace('(n+2)', '6').replace('(n-', '').replace('n ',
                                                                                                                  '4').replace(
                'n', '').replace('()', '').replace('(', '')
            try:
                metb_list = re.findall(r'(\d*\.?\d*[a-z]*)?\s*(\S+)', metb.replace(' ', ''))[0]
            except:
                pan=1
                break
            met_id = metb_list[1]
            try:
                metb_name = find_tarname('NAME', met_id)
                metb_weather = find_tarname2('FORMULA', met_id)[1]
                metb_formula = find_tarname2('FORMULA', met_id)[0]
                metb_ANNO = find_tarname('CHARGE', met_id)
                metbolic = cobra.Metabolite(
                    met_id,
                    name=metb_name,
                    formula=metb_formula,
                    compartment='c',
                    charge=metb_ANNO
                )
                try:
                    rule = model.metabolites[model_reaction.index(metb_name)].id
                    if model.metabolites[
                        model_reaction.index(metb_name)].compartment=='c':
                        try:
                            rule = model.metabolites[model_reaction.index(metb_name.lower())].id
                            if not metb_list[0]:
                                reac_metabolic[rule] = -1
                            elif metb_list[0] == 'n':
                                reac_metabolic[rule] = -4
                            else:
                                reac_metabolic[rule] = -eval(metb_list[0])
                        except:
                            model.add_metabolites([metbolic])
                            rule = model.metabolites.get_by_id(met_id).id
                            if not metb_list[0]:
                                reac_metabolic[rule] = -1
                            elif metb_list[0] == 'n':
                                reac_metabolic[rule] = -4
                            else:
                                reac_metabolic[rule] = -eval(metb_list[0])
                    else:
                        print(rule)
                        if not metb_list[0]:
                            reac_metabolic[rule] = -1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = -4
                        else:
                            reac_metabolic[rule] = -eval(metb_list[0])
                except:
                    try:
                        rule = model.metabolites[model_reaction.index(metb_name.lower())].id
                        if not metb_list[0]:
                            reac_metabolic[rule] = -1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = -4
                        else:
                            reac_metabolic[rule] = -eval(metb_list[0])
                    except:
                        model.add_metabolites([metbolic])
                        rule = model.metabolites.get_by_id(met_id).id
                        if not metb_list[0]:
                            reac_metabolic[rule] = -1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = -4
                        else:
                            reac_metabolic[rule] = -eval(metb_list[0])
            except Exception as e:
                pan = 1
                print('here is a wrong for {} {} {}'.format(met_id, metb, row['REACTION']))
                print(e)
                break
        for metb in metb_right:
            model_reaction = []
            for i in model.metabolites:
                model_reaction.append(i.name.lower())
            metb = metb.replace('(n-2)', '2').replace('(m)', '').replace('(n+2)', '6').replace('(n-', '').replace('n ',
                                                                                                                  '4').replace(
                'n', '').replace('()', '').replace('(', '')
            try:
                metb_list = re.findall(r'(\d*\.?\d*[a-z]*)?\s*(\S+)', metb.replace(' ', ''))[0]
            except:
                pan=1
                break
            met_id = metb_list[1]
            try:
                metb_name = find_tarname('NAME', met_id)
                metb_weather = find_tarname2('FORMULA', met_id)[1]
                metb_formula = find_tarname2('FORMULA', met_id)[0]
                metb_ANNO = find_tarname('CHARGE', met_id)
                metbolic = cobra.Metabolite(
                    met_id,
                    name=metb_name,
                    formula=metb_formula,
                    compartment='c',
                    charge=metb_ANNO
                )
                try:
                    rule = model.metabolites[model_reaction.index(metb_name)].id
                    if metb_formula != model.metabolites[
                        model_reaction.index(metb_name)].compartment=='c':
                        try:
                            rule = model.metabolites[model_reaction.index(metb_name.lower())].id
                            print(rule)
                            if not metb_list[0]:
                                reac_metabolic[rule] = 1
                            elif metb_list[0] == 'n':
                                reac_metabolic[rule] = 4
                            else:
                                reac_metabolic[rule] = eval(metb_list[0])
                        except:
                            model.add_metabolites([metbolic])
                            rule = model.metabolites.get_by_id(met_id).id
                            if not metb_list[0]:
                                reac_metabolic[rule] = 1
                            elif metb_list[0] == 'n':
                                reac_metabolic[rule] = 4
                            else:
                                reac_metabolic[rule] = eval(metb_list[0])
                    else:
                        print(rule)
                        if not metb_list[0]:
                            reac_metabolic[rule] = 1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = 4
                        else:
                            reac_metabolic[rule] = eval(metb_list[0])
                except:
                    try:
                        rule = model.metabolites[model_reaction.index(metb_name.lower())].id
                        print(rule)
                        if not metb_list[0]:
                            reac_metabolic[rule] = 1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = 4
                        else:
                            reac_metabolic[rule] = eval(metb_list[0])
                    except:
                        model.add_metabolites([metbolic])
                        rule = model.metabolites.get_by_id(met_id).id
                        if not metb_list[0]:
                            reac_metabolic[rule] = 1
                        elif metb_list[0] == 'n':
                            reac_metabolic[rule] = 4
                        else:
                            reac_metabolic[rule] = eval(metb_list[0])
            except Exception as e:
                pan = 1
                print('here is a wrong for {} {} {}'.format(met_id, metb, row['REACTION']))
                print(e)
                break
        if pan == 0:
            if check_reactions(reac_metabolic,model)==0:
                continue
            model.add_reactions([Reaction])
            Reaction.add_metabolites(reac_metabolic)
            Reaction.gene_reaction_rule = gpr.iat[list(gpr['reaction']).index(Reaction.id), 2]
    for i in model.reactions:
        values = np.array(list(i.metabolites.values()))
        if np.sum(values >= 0) == 0 and i.lower_bound != 0 and i.id not in model.medium:
            print(i.id)
            model.reactions.get_by_id(i.id).knock_out()
    cobra.io.save_yaml_model(model, f'models/tarmodel_{name}.yml')




