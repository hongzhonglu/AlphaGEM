import os
import pandas as pd
taryeast=pd.read_excel('ziyuan/Kpneumoniae.xlsx')
name='Kpneumoniae'
os.system('chmod 660 /home/pickleopear/.wget-hsts')
for i in range(len(taryeast.index)):
    if type(taryeast.iat[i,2])==float:
        continue
    cmd1=taryeast.iat[i,0]
    try:
        os.system(f'wget https://alphafold.ebi.ac.uk/files/AF-{cmd1}-F1-model_v4.pdb -P /home/pickleopear/PycharmProjects/Alpha-GEMs/struct_data/taryeast/{name} ')
    except:
        continue