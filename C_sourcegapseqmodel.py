import cobra
model = cobra.io.read_sbml_model(f'models/KP.xml')
medium = model.medium
medium['EX_cpd00023_e0']=0
medium['EX_cpd00027_e0']=0
medium['EX_cpd00041_e0']=0
medium['EX_cpd00382_e0']=0
medium['EX_cpd90020_e0']=0
medium['EX_cpd01133_e0']=0
medium['EX_cpd00059_e0']=0
medium['EX_cpd00220_e0']=0
medium['EX_cpd00311_e0']=0
medium['EX_cpd00098_e0']=0
medium['EX_cpd00305_e0']=0
import pandas as pd
table = pd.read_excel('ziyuan/TableS2.xls')
medium2 = medium
reac = []
carbons = pd.DataFrame()
for i in model.reactions:
    reac.append(i.name.lower())
for index, row in table.iterrows():
    medium2 = medium
    growth = pd.DataFrame({'Carbon': row['C source'], 'fba': 0, 'edata': row['data']},
                          index=[0])
    try:
        try:
            medium2[model.reactions[reac.index(row['C source'].lower() + '-e0 exchange')].id] = 1.0
        except:
            try:
                medium2[model.reactions[reac.index('d-' + row['C source'].lower() + '-e0 exchange')].id] = 1.0
            except:
                medium2[model.reactions[reac.index('l-' + row['C source'][2:].lower() + '-e0 exchange')].id] = 1.0
                print('special' + row['C source'][2:])
    except:
        print('not find ' + row['C source'])
        carbons = pd.concat([carbons, growth])
        continue
    model.medium = medium2
    growth['fba'] = model.optimize().objective_value
    carbons = pd.concat([carbons, growth])
    print(model.optimize().objective_value)
    try:
        try:
            medium2[model.reactions[reac.index(row['C source'].lower() + '-e0 exchange')].id] = 0
        except:
            try:
                medium2[model.reactions[reac.index('d-' + row['C source'].lower() + '-e0 exchange')].id] = 0
            except:
                medium2[model.reactions[reac.index('l-' + row['C source'].lower() + '-e0 exchange')].id] = 0
    except:
        continue
ymodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
for i in range(len(carbons.index)):
    if carbons.iat[i, 1] >= 0.001:
        carbons.iat[i, 1] = 'G'
    else:
        carbons.iat[i, 1] = 'NG'
carbons.to_excel('ccarbon2.xlsx')
carbons = pd.read_excel('ccarbon2.xlsx')
tp = 0;tn = 0;fn = 0;fp = 0
for i in range(len(carbons.index)):
    if carbons.iat[i, 2] == 'G' and carbons.iat[i, 3] == 'E-G':
        tp += 1
    if carbons.iat[i, 2] == 'NG' and carbons.iat[i, 3] == 'E-G':
        fn += 1
    if carbons.iat[i, 2] == 'NG' and carbons.iat[i, 3] == 'E-NG':
        tn += 1
    if carbons.iat[i, 2] == 'G' and carbons.iat[i, 3] == 'E-NG':
        fp += 1