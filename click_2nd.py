import pandas as pd
import numpy as np
juzhen2=pd.DataFrame()
taryeast=pd.read_excel('data_available/candida.xlsx')
juzhen1=pd.read_excel('juzhen/juzhen_homolog.xlsx')
kkkk=[]
for i in range(len(juzhen1.index)):
    kkkk.append(juzhen1.iat[i,2])
for i in range(len(taryeast.index)):
    pan=0
    print(i)
    if taryeast.iat[i,0] in kkkk:
            pan=1
    if pan==0:
            juzhen2=pd.concat([juzhen2,pd.DataFrame([taryeast.iat[i,0]])])
juzhen2.index=range(len(juzhen2.index))
juzhen2.to_excel('juzhen/juzhen_2ndkout.xlsx')