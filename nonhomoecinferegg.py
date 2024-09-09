import pandas as pd
def nonhome(name):
    kegg2rhea = pd.read_csv('ziyuan/rhea2kegg_reaction.tsv', sep='\t')
    keggdict = kegg2rhea.set_index('ID')['RHEA_ID'].to_dict()
    egg = pd.read_excel(f'juzhen/eggec2{name}.xlsx')
    eggec = pd.DataFrame()
    for genes in range(len(egg.index)):
        if egg.iat[genes, 3] == '-':
            continue
        for ec in egg.iat[genes, 3].split(','):
            try:
                eggec = pd.concat([eggec, pd.DataFrame({'reaction': egg.iat[genes, 1], 'ec': keggdict[ec]}, index=[1])])
            except:
                continue
    result = eggec.groupby('ec')['reaction'].apply(lambda x: ' or '.join(x)).reset_index()
    result.to_excel(f'juzhen/{name}gpregg.xlsx')


