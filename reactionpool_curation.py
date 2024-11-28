import cobra
import pandas as pd
from tqdm import tqdm
reactionpool=cobra.io.load_json_model('models/reaction_pool.json')
metnames=pd.read_excel('data_available/metanetxrefcpd.xlsx')
rnxnames=pd.read_excel('data_available/metanetxrefrxn.xlsx')
chebidict={}
for index,row in metnames.iterrows():
    for id in row['chebi'][1:-1].split(', '):
        chebidict[id[1:-1]]=row[0]
biggdict={row[0]:[] for index,row in metnames.iterrows()}
for index,row in metnames.iterrows():
    for id in row['bigg.metabolite'][1:-1].split(', '):
        biggdict[row[0]].append(id[1:-1])
keggdict={row[0]:[] for index,row in metnames.iterrows()}
for index,row in metnames.iterrows():
    for id in row['kegg.compound'][1:-1].split(', '):
        keggdict[row[0]].append(id[1:-1])


for met in tqdm(list(reactionpool.metabolites)):
    chebi=met.annotation['chebi']
    try:
        met.id="chebi"+chebi+'_c'
    except:
        reactions=list(met.reactions)
        for i in reactions:
            metlist=reactionpool.reactions.get_by_id(i.id).metabolites
            try:
                 metlist[reactionpool.metabolites.get_by_id("chebi"+chebi+'_c')]=metlist[reactionpool.metabolites.get_by_id("chebi"+chebi+'_c')]+metlist.pop(met)
            except:
                metlist[reactionpool.metabolites.get_by_id("chebi" + chebi + '_c')]=metlist.pop(met)
            reactionchange=cobra.Reaction(lower_bound=i.lower_bound, upper_bound=i.upper_bound,id=i.id,name=i.name)
            reactionchange.annotation=i.annotation
            reactionchange.add_metabolites(metlist)
            reactionpool.remove_reactions([i])
            reactionpool.add_reactions([reactionchange])
        reactionpool.remove_metabolites(met)
    try:
        met.annotation['bigg.metabolite']=biggdict[chebidict[chebi]]
    except:
        pass
    try:
        met.annotation['kegg.compound']=keggdict[chebidict[chebi]]
    except:
        pass
    try:
        met.annotation['metanetx'] = chebidict[chebi]
    except:
        pass


cobra.io.save_json_model(reactionpool,'models/reaction_pool_met.json')
reactionpool=cobra.io.load_json_model('models/reaction_pool_met.json')
for ex in list(reactionpool.exchanges):
    metabolic2=ex.products[0].copy()
    metabolic2.compartment='e'
    metabolic2.id=metabolic2.id.replace('_c','_e')
    ex.add_metabolites({metabolic2:-1})
cobra.io.save_json_model(reactionpool,'models/reaction_pool_met_trans.json')

reactionpool=cobra.io.load_json_model('models/reaction_pool_met_trans.json')
rheadict={}
for index,row in rnxnames.iterrows():
    for id in row['rheaR'][1:-1].split(', '):
        rheadict[id[1:-1]]=row[0]
keggdict={row[0]:[] for index,row in rnxnames.iterrows()}
for index,row in rnxnames.iterrows():
    for id in row['keggR'][1:-1].split(', '):
        keggdict[row[0]].append(id[1:-1])
biggdict={row[0]:[] for index,row in rnxnames.iterrows()}
for index,row in rnxnames.iterrows():
    for id in row['bigg.reaction'][1:-1].split(', '):
        biggdict[row[0]].append(id[1:-1])


for rea in tqdm(list(reactionpool.reactions)):
    rhea = rea.id.split('_')[1]
    try:
        rea.id=rea.id+'_c'
    except:
        print(rea.id)
    try:
        rea.annotation['bigg.reaction']=biggdict[rheadict[rhea]]
    except:
        pass
    try:
        rea.annotation['keggR']=keggdict[rheadict[rhea]]
    except:
        pass
    try:
        rea.annotation['metanetx']=rheadict[rhea]
    except:
        pass

cobra.io.save_json_model(reactionpool,'models/reaction_pool_final.json')
reactionpool=cobra.io.load_json_model('models/reaction_pool_final.json')
for met in reactionpool.reactions:
    met.id = met.id[:-2]
for met in reactionpool.metabolites:
    met.id=met.id.replace('_c','')
    met.compartment='part1'
for met in reactionpool.metabolites:
    if met.id[-2:]=='_e':
        met.compartment='part2'
for met in reactionpool.metabolites:
    if met.compartment=='part1':
        met.compartment='in'
    else:
        met.compartment='out'
for met in reactionpool.metabolites:
    if met.id[-2:] == '_e':
        met.id=met.id[:-2]+'_out'
    else:
        met.id=met.id+'_in'
cobra.io.save_json_model(reactionpool,'models/reaction_pool_noncompartment.json')



