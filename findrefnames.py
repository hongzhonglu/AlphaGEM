import pandas as pd
from config import refname
df=pd.read_excel(f'ziyuan/{refname}.xlsx')
entrys=[]
for index,row in df.iterrows():
    entrys.append(row.iloc[0])
modelnames=[]
for index,row in df.iterrows():
    modelnames.append(row.iloc[2])
def findmodelname(modelname):
    try:
        return df.iat[modelnames.index(modelname),0]
    except:
        print(f'This gene({modelname}) is not in the table ,please give me the right table')
def findentryname(entry):
    try:
        return df.iat[entrys.index(entry),2]
    except:
        print(f'This gene({entry}) is not in the table ,please give me the right table')
