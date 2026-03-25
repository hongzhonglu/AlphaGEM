import pandas as pd
import cobra
import pickle
from multiple_annotition_for_nonhomogene import calculate_score

def filter_plmsearch(merged,general2spe,filter_num,done=True):
    if done:
        results = []
        for index, row in merged.iterrows():
            if str(row['rhea']) in general2spe:
                for i in general2spe[str(row['rhea'])]:
                    safe=0
                    if len(merged[(merged['reaction'] == row['reaction']) & (merged['rhea'] == int(i))])==0:
                        continue
                    else:
                        for _,m in merged[(merged['reaction'] == row['reaction']) & (merged['rhea'] == int(i))].iterrows():
                            if m['anno_D']=='plmsearch':
                                safe=1
                                continue
                    if safe==0:
                        continue
                    row['rhea'] = i
                    results.append(row)
            else:
                    results.append(row)
        merged = pd.concat(results, axis=1).T.drop_duplicates()
        merged = merged.groupby(['reaction', 'rhea'], as_index=False).first()
        merged['final_score'] = merged.apply(calculate_score, axis=1)
        resultec=merged[merged['final_score'] >= filter_num].reset_index(drop=True)
        return resultec
    else:
        results=[]
        for index,row in merged.iterrows():
            if str(row['rhea']) in general2spe:
                for i in general2spe[str(row['rhea'])]:
                    row['rhea'] = i
                    results.append(row)
            else:
                results.append(row)
        merged = pd.concat(results,axis=1).T.drop_duplicates()
        merged = merged.groupby(['reaction', 'rhea'], as_index=False).first()
        merged['final_score'] = merged.apply(calculate_score, axis=1)
        resultec = merged[merged['final_score'] >= filter_num].reset_index(drop=True)
        return resultec

def get_spe_reactions(name,filter_num=0.4):
    merged=pd.read_excel(f'working/{name}/{name}multiple_score.xlsx').drop_duplicates()
    pklopen = open('data_available/general2specific2.pkl', 'rb')
    general2spe = pickle.load(pklopen)
    resultec=filter_plmsearch(merged,general2spe,filter_num,done=True)
    result = resultec.groupby('rhea')['reaction'].apply(lambda x: ' or '.join(x)).reset_index()
    result.to_excel(f'working/{name}/{name}gpr_score.xlsx')