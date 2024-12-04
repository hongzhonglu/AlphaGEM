import numpy as np
import pandas as pd
import tqdm
def dropblank(key):
    if len(key[6])+len(key[10])+len(key[4])==0:
        return np.nan
    else:
        return key[0]


def dropblank2(key):
    if len(key[0])+len(key[6])+len(key[12])==0:
        return np.nan
    else:
        return key[0]


data=pd.read_csv('data_available/chem_xref.tsv', sep='\t')
dict={}
for index,row in tqdm.tqdm(data.iterrows()):
    try:
        dict[row['ID']][row['source'].split(':')[0]].append(row['source'].split(':')[1])
    except:
        dict[row['ID']]={'ID':[], 'mnx':[], 'seed.compound':[], 'seedM':[], 'bigg.metabolite':[], 'biggM':[],
       'chebi':[], 'envipath':[], 'envipathM':[], 'hmdb':[], 'kegg.compound':[], 'keggC':[],
       'metacyc.compound':[], 'metacycM':[], 'reactome':[], 'reactomeM':[],
       'sabiork.compound':[], 'sabiorkM':[], 'CHEBI':[], 'SLM':[], 'lipidmaps':[],
       'lipidmapsM':[], 'slm':[], 'keggD':[], 'kegg.drug':[], 'kegg.glycan':[], 'keggG':[],
       'rheaG':[]}
    dict[row['ID']]['ID']=row['ID']
df=pd.DataFrame(dict).T

df['ID']=df.apply(dropblank,axis=1)
# df2=df.drop(list(df.apply(dropblank,axis=1)))
df2=df.dropna(axis=0,subset=['ID'])
df2.to_excel('data_available/metanetxrefcpdALL.xlsx')
# step0，将rhea反应mapping到father反应中
# step1，以matnetx为标准进行xref标准化
# step2，将reference进行MANX标准化
# step3，检查反应是否存在，若存在，加入基因关系用model——build的函数
# step4，检查代谢物存不存在，若不存在，添加反应后，获取反应中代谢物dict，删除dict，加入模型，重构dict，完成任务？（待定实现方式，或许可以统一化id，若没有再加入代谢物）
data=pd.read_csv('data_available/reac_xref.tsv', sep='\t')
dict={}
for index,row in tqdm.tqdm(data.iterrows()):
    try:
        dict[row['ID']][row['source'].split(':')[0]].append(row['source'].split(':')[1])
    except:
        dict[row['ID']]={'bigg.reaction':[],'biggR':[],'metacyc.reaction':[],'metacycR':[],'mnx':[],'rhea':[],'rheaR':[],'seed.reaction':[],'seedR':[],'sabiork.reaction':[],'sabiorkR':[],'kegg.reaction':[],'keggR':[]}
df=pd.DataFrame(dict).T
df['ID']=df.apply(dropblank2,axis=1)
df2=df.dropna(axis=0,subset=['ID'])
df2.to_excel('data_available/metanetxrefrxnALL.xlsx')