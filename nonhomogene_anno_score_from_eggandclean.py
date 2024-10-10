import pandas as pd
from rdflib import Graph, Namespace, URIRef
import cobra


def calculate_score(row):
    if pd.notnull(row['anno_B']):
        return row['score'] * 0.85 + 0.2
    else:
        return row['score'] * 0.85


def findtargetreaction(g, ec):
    rhea=[]
    rh = Namespace("http://rdf.rhea-db.org/")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    reaction_query = f"""
       SELECT DISTINCT ?x ?reaction
       WHERE {{
         ?x rdfs:subClassOf rh:Reaction ;
            rh:ec <http://purl.uniprot.org/enzyme/{ec}> ;
            rh:bidirectionalReaction ?y ;
       }}
       LIMIT 100
       """
    reactions = g.query(reaction_query)
    rhea=[str(row.x).split("/")[-1] for row in reactions]
    return rhea

def nonhome(name,clean_use):
    g=Graph()
    g.parse('rhea.rdf', format='xml')
    kegg2rhea = pd.read_csv('ziyuan/rhea2kegg_reaction.tsv', sep='\t')
    keggdict = kegg2rhea.set_index('ID')['MASTER_ID'].to_dict()
    egg = pd.read_excel(f'juzhen/eggec2{name}.xlsx')
    eggec = pd.DataFrame()
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
    cleanecrhea = pd.DataFrame()
    if clean_use:
        cleanecs = pd.DataFrame()
        clean = pd.read_csv(f'CLEAN/app/results/inputs/{name}_homoleft_maxsep.csv', sep='\t', header=None,
                            names=['clean'])
        for i in range(len(clean.index)):
            cleanecs = pd.concat([cleanecs, pd.DataFrame(
                {'reaction': clean.iat[i, 0].split('|')[1], 'ec': clean.iat[i, 0][clean.iat[i, 0].find('EC:'):]},
                index=[1])])
        for genes in range(len(cleanecs.index)):
            if cleanecs.iat[genes, 1] == '-':
                continue
            for ec in cleanecs.iat[genes, 1].split(','):
                cleanec = pd.concat([cleanec,
                                     pd.DataFrame({'reaction': cleanecs.iat[genes, 0], 'ec': ec.split('/')[0][3:],'score':float(ec.split('/')[1])},
                                                  index=[1])])
        for i in range(len(cleanec.index)):
            for rheas in findtargetreaction(g,cleanec.iat[i,1]):
                cleanecrhea = pd.concat([cleanecrhea,
                                     pd.DataFrame({'reaction': cleanec.iat[i, 0], 'rhea':str(rheas) ,'anno':'clean','score':cleanec.iat[i, 2]},
                                                  index=[1])])
        merged = pd.merge(cleanecrhea, eggec, how='outer', on=['reaction', 'rhea'], suffixes=('_A', '_B'))
        merged['final_score'] = merged.apply(calculate_score, axis=1)
        merged['final_score'].fillna(0.75, inplace=True)
        reaction_counts = merged['reaction'].value_counts()
        merged['final_score'] += merged['reaction'].apply(lambda x: 0.3 if reaction_counts[x] == 1 else 0)
        merged.to_excel(f'juzhen/{name}eggnog&clean_score.xlsx', index=False)
        resultec=merged[merged['final_score'] >= 0.75].reset_index(drop=True)
    result = resultec.groupby('rhea')['reaction'].apply(lambda x: ' or '.join(x)).reset_index()
    result.to_excel(f'juzhen/{name}gpr_score.xlsx')
