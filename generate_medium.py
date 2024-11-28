import cobra
import pandas as pd
from fuzzywuzzy import fuzz,process
import pickle

model=cobra.io.read_sbml_model('models/yeast-GEM.xml')
medium=pd.read_excel('data_available/mediumfull.xls')
exchanges={i.name.split(' exchange')[0].lower():i.id for i in model.exchanges}
modelmedium=model.medium
modelmedium.pop('r_1714')
mediumleft={}

for index, row in medium.iterrows():
    carbon,score=process.extractOne(row['name'].lower(),exchanges.keys())
    if row['formula']==model.reactions.get_by_id(exchanges[carbon]).reactants[0].formula:
        print(carbon,score,row['name'])
        modelmedium[exchanges[carbon]]=5.0
    else:
        booll=False
        for exchang in exchanges.values():
            if row['formula'] == model.reactions.get_by_id(exchang).reactants[0].formula:
                booll=True
                mediumleft[row['name']]=[carbon,score,row['formula'],model.reactions.get_by_id(exchang).reactants[0].name]
                medium[exchang]=5.0
        if not booll:
            mediumleft[row['name']]=[carbon,score,exchanges[carbon]]
    #print(carbon,score,row['name'])
modelmedium['r_2028']=5.0
mediumpkl=open('data_available/yeast_full_medium.pkl','wb')
pickle.dump(modelmedium,mediumpkl)
mediumpkl.close()
mediumpkl=open('data_available/yeast_full_medium.pkl','rb')
mediumpk=pickle.load(mediumpkl)
for reaction in modelmedium.keys():
    print(model.reactions.get_by_id(reaction).name)

model=cobra.io.load_json_model('models/iML1515.json')
medium=pd.read_excel('data_available/mediumfull.xls')
exchanges={i.name.split(' exchange')[0].lower():i.id for i in model.exchanges}
modelmedium=model.medium
modelmedium.pop('EX_glc__D_e')
mediumleft={}
for m in modelmedium.keys():
    print(model.reactions.get_by_id(m).reactants[0].name,model.reactions.get_by_id(m).reactants[0].formula)
for index, row in medium.iterrows():
    carbon,score=process.extractOne(row['name'].lower(),exchanges.keys())
    if row['formula']==model.reactions.get_by_id(exchanges[carbon]).reactants[0].formula:
        modelmedium[exchanges[carbon]]=5.0
    else:
        booll=False
        for exchang in exchanges.values():
            if row['formula'] == model.reactions.get_by_id(exchang).reactants[0].formula:
                print(model.reactions.get_by_id(exchang).reactants[0].name,row['name'])
                booll=True
                mediumleft[row['name']]=[carbon,score,row['formula'],model.reactions.get_by_id(exchang).reactants[0].name]
                medium[exchang]=5.0
        if not booll:
            mediumleft[row['name']]=[carbon,score,exchanges[carbon]]

    #print(carbon,score,row['name'])
metafl=model.metabolites.ribflv_c.copy()
metafl.id='ribflv_e'
metafl.compartment='e'
Reac=cobra.Reaction('EX_ribflv_e')
Reac.add_metabolites({metafl:-1})
Reac.bounds=[0,1000]
model.add_reactions([Reac])
cobra.io.save_json_model(model,'models/iML1515.json')
modelmedium['EX_ribflv_e']=5.0
mediumpkl=open('data_available/ecoli_full_medium.pkl','wb')
pickle.dump(modelmedium,mediumpkl)
mediumpkl.close()
mediumpkl=open('data_available/ecoli_full_medium.pkl','rb')
mediumpk=pickle.load(mediumpkl)
for reaction in modelmedium.keys():
    print(model.reactions.get_by_id(reaction).name)