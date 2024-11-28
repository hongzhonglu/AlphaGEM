import cobra
import pandas as pd
modelyeast=cobra.io.read_sbml_model('models/yeast-GEM-old.xml')
metnames=pd.read_excel('data_available/metanetxrefcpd.xlsx')
rnxnames=pd.read_excel('data_available/metanetxrefrxn.xlsx')
reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
chebi2metdict1={}
keggdict={}
for index,row in metnames.iterrows():
    for id in row['kegg.compound'][1:-1].split(', '):
        keggdict[id[1:-1]]=row['chebi'][1:-1].split(', ')[0][-1:1]
metdict={}
for index,row in metnames.iterrows():
        metdict[row[0]]=row['chebi'][1:-1].split(', ')[0][-1:1]
for me in reactionpool.metabolites:
    try:
       chebi2metdict1[me.annotation['metanetx']]=me.annotation['chebi']
    except:
        pass
chebi2keggdict={}
for me in reactionpool.metabolites:
    try:
        for kegg in me.annotation['kegg.compound']:
            chebi2keggdict[kegg]=me.annotation['chebi']
    except:
        pass
count=0
mets1=[]
for met in modelyeast.metabolites:
    try:
        met.id='chebi'+met.annotation['chebi'].split(':')[1]+'_'+met.compartment
    except:
        count+=1
        mets1.append(met)
modelyeast.metabolites.s_0450.id='biomass'
mets1.remove(modelyeast.metabolites.biomass)
count2=0
mets2=[]
for met in mets1:
    try:
        met.id='chebi'+chebi2metdict1[met.annotation['metanetx.chemical']]+'_'+met.compartment
    except:
        count2+=1
        mets2.append(met)
count3=0
mets3=[]
for met in mets2:
    try:
        met.id='chebi'+chebi2keggdict[met.annotation['kegg.compound']]+'_'+met.compartment
    except:
        print(met)
        count3+=1
        mets3.append(met)

rhea2metanetxdict={}
for reaction in reactionpool.reactions:
    try:
        rhea2metanetxdict[reaction.annotation['metanetx']]=reaction.id
    except:
        pass
rhea2keggdict={}
for reaction in reactionpool.reactions:
    try:
        for kegg in reaction.annotation['keggR']:
             rhea2keggdict[kegg]=reaction.id
    except:
        pass

reacts=[]
count=0
for r in modelyeast.reactions:
    try:
        if len(r.compartments) ==1:
           r.id=rhea2metanetxdict[r.annotation['metanetx.reaction']][:-2]+'_'+list(r.compartments)[0]
        else:
            count += 1
            reacts.append(r)
    except:
        try:
            if len(r.compartments) == 1:
                r.id = rhea2keggdict[r.annotation['kegg.reaction']][:-2]+'_'+list(r.compartments)[0]
            else:
                count += 1
                reacts.append(r)
        except:
            count += 1
            reacts.append(r)
cobra.io.write_sbml_model(modelyeast,'models/yeast-GEM.xml')


