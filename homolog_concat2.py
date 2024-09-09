import pandas as pd
def homo(name):
    gx1 = pd.read_excel(f'juzhen/juzhen_model1{name}.xlsx')
    gx4 = pd.read_excel(f'juzhen/juzhen_model3{name}.xlsx')
    gx3 = pd.DataFrame()
    for i in range(len(gx1.index)):
        gx3 = pd.concat([gx3, pd.DataFrame({
            0: [gx1.iat[i, 1]],
            1: [gx1.iat[i, 2]]
        })])
    for i in range(len(gx4.index)):
        try:
            a=gx3[0]
            if gx3.iat[list(a).index(gx4.iat[i,1]),1]==gx4.iat[i,2]:
                continue
        except:
            a=0
        gx3 = pd.concat([gx3, pd.DataFrame({
            0: [gx4.iat[i, 1]],
            1: [gx4.iat[i, 2]]
        })])
    gx3.index = range(len(gx3.index))
    gx3.to_excel(f'juzhen/juzhen_homolog{name}.xlsx')