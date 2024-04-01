import pandas as pd
import csv
import numpy as np
with open('ziyuan/chem_xref.tsv', 'r') as f:
     reader = csv.reader(f)
     row2=[]
     che=''
     for row in reader:
         print(row)
         if 'kegg.compound' in row[0]:
             row2.append(','.join(row))
rxn=pd.DataFrame()
for i in row2:
    rxn=pd.concat([rxn,np.transpose(pd.DataFrame(i.replace('kegg.compound:','').split('\t')))])
with open('ziyuan/chem_prop.tsv','r') as f:
    reader=csv.reader(f)
    row3=[]
    for row in reader:
        if 'keggC' in ','.join(row):
          row3.append(','.join(row).replace('keggC:',''))
    row31=pd.DataFrame()
    for row in row3:
        row31 = pd.concat([row31, np.transpose(pd.DataFrame(row.split('\t')))])
        print(row)
row31.drop(5,axis=1,inplace=True)
row31.drop(6,axis=1,inplace=True)
row31.drop(7,axis=1,inplace=True)
row31.drop(8,axis=1,inplace=True)
row31.index=range(len(row31.index))
for ind,row in row31.iterrows():
    if row[4]=='':
        print(ind)
        row31=row31.drop(ind)
row31.to_excel('ziyuan/cpd.xlsx')
row31.columns=['0','NAME','ID','FORMULA','CHARGE']