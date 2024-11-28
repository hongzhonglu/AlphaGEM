import pandas as pd
from rdflib import Graph, Namespace, URIRef
import cobra
import os


def calculate_score(row):
    score=0
    if pd.notnull(row['anno_B']):
        score+=0.2
    if pd.notnull(row['anno_A']):
        score+=0.2
    if pd.notnull(row['anno_C']):
        score+=0.2
    if pd.notnull(row['anno_D']):
        score+=0.2
    return score


def findtargetreaction(g, ec):
    rhea=[]
    rh = Namespace("http://rdf.rhea-db.org/")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    reaction_query = f"""
       SELECT DISTINCT ?x ?reaction
       WHERE {{
         ?x rdfs:subClassOf rh:Reaction ;
            rh:ec <http://purl.uniprot.org/enzyme/{ec}> ; ;
       }}
       LIMIT 100
       """
    reactions = g.query(reaction_query)
    rhea=[str(row.x).split("/")[-1] for row in reactions]
    return rhea

def nonhome(name,clean_use,deepec_use,plm_use=True):
    g=Graph()
    g.parse('rhea.rdf', format='xml')
    kegg2rhea = pd.read_csv('data_available/rhea2kegg_reaction.tsv', sep='\t')
    keggdict = kegg2rhea.set_index('ID')['MASTER_ID'].to_dict()
    egg = pd.read_excel(f'working/{name}/eggec2{name}.xlsx')
    eggec = pd.DataFrame(data=None,columns=['reaction','rhea','anno'])
    for genes in range(len(egg.index)):
        if egg.iat[genes, 3] == '-':
            continue
        for ec in egg.iat[genes, 3].split(','):
            try:
                eggec = pd.concat(
                    [eggec, pd.DataFrame({'reaction': egg.iat[genes, 1], 'rhea': str(keggdict[ec]),'anno':'eggnog_mapper'}, index=[1])])
            except:
                continue

    resultec=eggec
    cleanec = pd.DataFrame()
    cleanecrhea = pd.DataFrame(data=None,columns=['reaction','rhea','anno'])
    deepec = pd.DataFrame()
    deepecrhea = pd.DataFrame(data=None,columns=['reaction','rhea','anno'])
    plmrhea = pd.DataFrame(data=None,columns=['reaction','rhea','anno'])
    filter_num=0.2
    if clean_use:
        filter_num=0.4
        cleanecs = pd.DataFrame()
        clean = pd.read_csv(f'working/{name}/{name}_homoleft_maxsep.csv', sep='\t', header=None,
                            names=['clean'])
        for i in range(len(clean.index)):
            cleanecs = pd.concat([cleanecs, pd.DataFrame(
                {'reaction': clean.iat[i, 0].split('|')[1], 'ec': clean.iat[i, 0][clean.iat[i, 0].find('EC:'):]},
                index=[1])])
        for genes in range(len(cleanecs.index)):
            if cleanecs.iat[genes, 1] == '-':
                continue
            for ec in cleanecs.iat[genes, 1].split(','):
                if float(ec.split('/')[1]) <=0.7:
                    continue
                cleanec = pd.concat([cleanec,
                                     pd.DataFrame({'reaction': cleanecs.iat[genes, 0], 'ec': ec.split('/')[0][3:]},
                                                  index=[1])])
        for i in range(len(cleanec.index)):
            for rheas in findtargetreaction(g,cleanec.iat[i,1]):
                cleanecrhea = pd.concat([cleanecrhea,
                                     pd.DataFrame({'reaction': cleanec.iat[i, 0], 'rhea':str(rheas) ,'anno':'clean'},
                                                  index=[1])])
    if deepec_use:
        filter_num=0.4
        path=os.path.join(os.getcwd(),'working',f'{name}',f'{name}_deepec_result','DeepECv2_result.txt')
        deepec=pd.read_csv(path,sep='\t')
        deepec['sequence_ID']=deepec['sequence_ID'].apply(lambda x: x.split('|')[1])
        deepec.fillna('',inplace=True)
        for genes in range(len(deepec.index)):
            if deepec.iat[genes, 1]=='':
                continue
            for ec in str(deepec.iat[genes, 1]).split(';'):
                if ec.endswith('-'):
                    continue
                for rheas in findtargetreaction(g,ec[3:]):
                    deepecrhea = pd.concat([deepecrhea,
                                     pd.DataFrame({'reaction': deepec.iat[genes, 0], 'rhea': str(rheas),'anno':'deepec'},
                                                  index=[1])])
    if plm_use:
        filter_num = 0.4
        plmresult=pd.read_excel(f'working/{name}/{name}_rheaid_plmsearch.xlsx')
        plmrhea['reaction']=plmresult['Query'].apply(lambda x:x.split('|')[1])
        plmrhea['rhea']=plmresult['Rhea_id']
        plmrhea['anno']='plmsearch'
        plmrhea=plmrhea.drop_duplicates()
        plmrhea['rhea']=plmrhea['rhea'].apply(lambda x:str(x))
        #for genes in range(len(plmresult.index)):
            #plmrhea = pd.concat([plmrhea,pd.DataFrame({'reaction':plmresult.iat[genes, 1], 'rhea':plmresult.iat[genes,2] ,'anno':'plm'},index=[1])])
    merged = pd.merge(cleanecrhea, eggec, how='outer', on=['reaction', 'rhea'], suffixes=('_A', '_B'))
    deepecrhea = deepecrhea.add_suffix('_C')
    deepecrhea.rename(columns={'reaction_C': 'reaction', 'rhea_C': 'rhea'}, inplace=True)
    merged = pd.merge(merged, deepecrhea, how='outer', on=['reaction', 'rhea'])
    plmrhea=plmrhea.add_suffix('_D')
    plmrhea.rename(columns={'reaction_D': 'reaction','rhea_D':'rhea'}, inplace=True)
    merged=pd.merge(merged, plmrhea, how='outer', on=['reaction', 'rhea'])
    merged['final_score'] = merged.apply(calculate_score, axis=1)
    resultec = merged[merged['final_score'] >= filter_num].reset_index(drop=True)
    merged.to_excel(f'working/{name}/{name}multiple_score.xlsx', index=False)
    result = resultec.groupby('rhea')['reaction'].apply(lambda x: ' or '.join(x)).reset_index()
    result.to_excel(f'working/{name}/{name}gpr_score.xlsx')