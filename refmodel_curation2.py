import cobra
import pandas as pd
metnames=pd.read_excel('data_available/metanetxrefcpd.xlsx')
rnxnames=pd.read_excel('data_available/metanetxrefrxn.xlsx')
reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
modelecoli=cobra.io.load_json_model('models/iML1515_old.json')

chebi2biggdict={}
for met in reactionpool.metabolites:
    try:
       for i in met.annotation['bigg.metabolite']:
           chebi2biggdict[i]=met.annotation['chebi']
    except:
        pass

count=0
for me in modelecoli.metabolites:
    id=me.id[0:-2]
    try:
        me.annotation['bigg.metabolite'] = id
        if me.elements==reactionpool.metabolites.get_by_id('chebi' + chebi2biggdict[id] + '_c').elements:
            me.id = 'chebi' + chebi2biggdict[id] + '_' + me.compartment
            me.annotation['chebi'] = chebi2biggdict[id]
        else:
            count += 1
    except KeyError as e:
        count+=1

rhea2biggdict={}
for rxn in reactionpool.reactions:
    try:
     for bigg in rxn.annotation['bigg.reaction']:
        rhea2biggdict[bigg]=rxn.id[:-2]
    except:
        pass
count=0
for reactions in modelecoli.reactions:
    reactions.annotation['bigg.reaction'] = reactions.id
    try:
        id = reactions.id
        reactions.annotation['rheaR'] = rhea2biggdict[id]
        if(len(reactions.compartments)==1):
            reactions.id=rhea2biggdict[id] + '_' + list(reactions.compartments)[0]
        else:
            count+=1
    except KeyError as e:
        count+=1
cobra.io.save_json_model(modelecoli, 'models/iML1515.json')