import os
import pandas as pd
import shutil
df=pd.read_excel('data_available/human.xlsx')
for index,row in df.iterrows():
    try:
       shutil.move(f'struct_data/human/AF-{row[0]}-F1-model_v4.pdb',f'struct_data/human2/AF-{row[0]}-F1-model_v4.pdb')
    except:
       continue