import os
import pandas as pd
taryeast=pd.read_excel('data_available/kpastoris.xlsx')
name='kpastoris'
#os.system('chmod 660 /home/pickleopear/.wget-hsts')
for i in range(len(taryeast.index)):
    if type(taryeast.iat[i,2])==float:
        continue
    cmd1=taryeast.iat[i,0]
    if os.path.exists(f'/dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/taryeast/{name}/AF-{cmd1}-F1-model_v4.pdb'):
        continue
    try:
        os.system(f'wget https://alphafold.ebi.ac.uk/files/AF-{cmd1}-F1-model_v4.pdb -P /dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/taryeast/{name}')
    except:
        continue
for path in os.listdir(f'/dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/taryeast/{name}'):
    if path.endswith('.1'):
        os.remove(f'/dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/taryeast/{name}/{path}')