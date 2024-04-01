import cobra
import pandas as pd
def modelbuild(refmodel,name):
    tarname = name
    tarmod = cobra.io.read_sbml_model(f'models/tarmod{name}.xml')
    if refmodel=='yeast-GEM.xml' or refmodel=='Sco-GEM.xml' or refmodel=='Human-GEM.xml':
      ymod = cobra.io.read_sbml_model(f'models/{refmodel}')
    if refmodel=='iML1515.json':
        ymod=cobra.io.load_json_model(f'models/{refmodel}')
    tarmodel = cobra.Model('Target_model')
    juzhen1 = pd.read_excel('juzhen/juzhen_homolog.xlsx')
    juzhen2 = []
    juzhen3 = []
    juzhen4 = []
    for i in range(len(ymod.reactions)):
        o = 0
        for j in range(len(tarmod.reactions)):
            if ymod.reactions[i].name == tarmod.reactions[j].name:
                o = 1
        if o == 0:
            juzhen2.append(i)
    for i in juzhen2:
        pp = 0
        for j in range(len(juzhen1.index)):
            try:
                if juzhen1.iat[j, 1] in ymod.reactions[i].gene_reaction_rule and juzhen1.iat[j,1]!=juzhen1.iat[j-1,1]:
                    pp += 1
            except:
                continue
        if pp >= 0.5 * len(ymod.reactions[i].genes) and len(ymod.reactions[i].genes)>=3:
            print(ymod.reactions[i].name)
            juzhen3.append(ymod.reactions[i].gene_reaction_rule)
            juzhen4.append(i)
            tarmod.add_reactions([ymod.reactions[i]])

    for i in range(len(juzhen4)):
        o = 0
        stra = ' and '
        for j in range(len(juzhen1.index)):
            try:
                if juzhen1.iat[j, 1] in juzhen3[i]:
                    o += 1
                    if o == 1:
                        tarmod.reactions.get_by_id(ymod.reactions[juzhen4[i]].id).gene_reaction_rule = juzhen1.iat[j, 2]
                    elif o > 1:
                        aa = [tarmod.reactions.get_by_id(ymod.reactions[juzhen4[i]].id).gene_reaction_rule,
                              juzhen1.iat[j, 2]]
                        tarmod.reactions.get_by_id(ymod.reactions[juzhen4[i]].id).gene_reaction_rule = stra.join(aa)
            except:
                continue

    for i in range(len(tarmod.reactions)):
        tarmodel.add_reactions([tarmod.reactions[i]])
    tarmodel.objective=ymod.objective
    fba = tarmodel.optimize()
    cobra.io.save_yaml_model(tarmodel, f'models/tarmodel{tarname}.yml')
