import pandas as pd
import cobra
import numpy as np
df=pd.read_excel('juzhen/whytarnsportor.xlsx')
model=cobra.io.read_sbml_model('models/yeast-GEM.xml')
genes=[]
df2=pd.DataFrame()
for i in model.genes:
    genes.append(i.id)
df.columns=['names','numbers','tmscore','yes/not']
for index,raw in df.iterrows():
    try:
        a=genes.index(raw['names'])
        df2=pd.concat([df2,np.transpose(pd.DataFrame(raw))])
    except:
        continue
df2=df2.sort_values('yes/not')
df2.to_excel('juzhen/trans.xlsx')