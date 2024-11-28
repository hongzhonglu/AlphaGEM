import cobra
import pandas as pd
metnames=pd.read_excel('data_available/metanetxrefcpd.xlsx')
rnxnames=pd.read_excel('data_available/metanetxrefrxn.xlsx')
reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
modelhuman=cobra.io.read_sbml_model('models/human-GEM-old.xml')
chebi2metdict={}
for met in reactionpool.metabolites:
    try:
       chebi2metdict[met.annotation['metanetx']]=met.annotation['chebi']
    except:
        pass

count=0
mets=[]
for met in modelhuman.metabolites:
    id=met.id
    try:
        met.id='chebi'+met.annotation['chebi'].split(':')[1]+'_'+met.compartment
    except:
        count+=1
        mets.append(met)


for meta in mets:
    id=meta.id
    try:
        meta.id='chebi'+chebi2metdict[meta.annotation['metanetx.chemical']]+'_'+meta.compartment
        meta.annotation['chebi']=chebi2metdict[id]
    except:
        count+=1

rhea2metdict={}
for rxn in reactionpool.reactions:
    try:
        rhea2metdict[rxn.annotation['metanetx']]=rxn.id[:-2]
    except:
        pass
rhea2keggdict={}
for rxn in reactionpool.reactions:
    try:
        for kegg in rxn.annotation['keggR']:
            rhea2keggdict[kegg]=rxn.id[:-2]
    except:
        pass
count=0
reacts2=[]
for reaction in modelhuman.reactions:
    try:
        if len(reaction.compartments)==1:
            reaction.id='rhea_'+reaction.annotation['rhea']+'_'+list(reaction.compartments)[0]
        else:
            reacts2.append(reaction)
            count+=1
    except:
        reacts2.append(reaction)
        count+=1

count=0
react=[]
for reac in reacts2:
    id=reac.id
    try:
        reac.annotation['rheaR']=rhea2metdict[reac.annotation['metanetx.reaction']]
        if len(reac.compartments) == 1:
            reac.id=rhea2metdict[reac.annotation['metanetx.reaction']]+'_'+list(reac.compartments)[0]
        else:
            count+=1
    except:
        count+=1
        react.append(reac)
count=0
for reacts in react:
    id=reacts.id
    try:
        reacts.annotation['rheaR']=rhea2keggdict[reacts.annotation['kegg.reaction']]
        if len(reacts.compartments) == 1:
            reacts.id=rhea2keggdict[reacts.annotation['kegg.reaction']]+'_'+list(reacts.compartments)[0]
        else:
            count+=1
    except:
        count+=1
cobra.io.write_sbml_model(modelhuman,'models/human-GEM.xml')
