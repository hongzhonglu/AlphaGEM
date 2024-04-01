import cobra
import pandas as pd
def C_source():
    name='candida'
    model = cobra.io.read_sbml_model(f'models/tarmodel__{name}text.xml')
    medium = model.medium
    medium['r_1714'] = 0
    import pandas as pd
    table = pd.read_excel('ziyuan/TableS1.xls')
    medium2 = medium
    reac = []
    for i in model.reactions:
        reac.append(i.name)
    for index, row in table.iterrows():
        medium2 = medium
        try:
            try:
                medium2[model.reactions[reac.index(row['C source'] + ' exchange')].id] = 1.0
            except:
                try:
                    medium2[model.reactions[reac.index(row['C source'].lower() + ' exchange')].id] = 1.0
                except:
                    medium2[model.reactions[reac.index(row['C source'][2:].lower() + ' exchange')].id] = 1.0
                    print('special' + row['C source'][2:])
        except:
            print('not find' + row['C source'])
            continue
        model.medium = medium2
        print(model.optimize().objective_value)
        try:
            try:
                medium2[model.reactions[reac.index(row['C source'] + ' exchange')].id] = 0
            except:
                try:
                    medium2[model.reactions[reac.index(row['C source'].lower() + ' exchange')].id] = 0
                except:
                    medium2[model.reactions[reac.index(row['C source'][2:].lower() + ' exchange')].id] = 0
        except:
            continue
    ymodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
C_source()
