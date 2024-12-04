import pandas as pd
def nonhome(name,clean_use):
    egg = pd.read_excel(f'juzhen/eggec2{name}.xlsx')
    eggec = pd.DataFrame()
    for genes in range(len(egg.index)):
        if egg.iat[genes, 2] == '-':
            continue
        for ec in egg.iat[genes, 2].split(','):
            eggec = pd.concat([eggec, pd.DataFrame({'reaction': egg.iat[genes, 1], 'ec': ec}, index=[1])])
    cleanec = pd.DataFrame()
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
                if float(ec.split('/')[1]) <= 0.8:
                    continue
                cleanec = pd.concat([cleanec,
                                     pd.DataFrame({'reaction': cleanecs.iat[genes, 0], 'ec': ec.split('/')[0][3:]},
                                                  index=[1])])
    resultec = pd.concat([cleanec, eggec], join='outer', ignore_index=True)
    result = resultec.groupby('ec')['reaction'].apply(lambda x: ' or '.join(x)).reset_index()
    result.to_excel(f'juzhen/{name}gpregg_clean.xlsx')


