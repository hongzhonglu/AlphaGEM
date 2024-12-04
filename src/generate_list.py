import pandas as pd
from Bio import SeqIO
import os
import shutil

def generare_list_without_structure(name):
    list=pd.DataFrame()
    genes=f'working/{name}/{name}.fasta'
    records = SeqIO.parse(genes, 'fasta')
    for record in records:
        list=pd.concat([list,pd.DataFrame({'Entry':[record.id.split('|')[1]],'Entry Name':[record.id.split('|')[2]],'Structure Data':[record.id.split('|')[1]]})])
    list.to_excel(f'working/{name}/{name}.xlsx', index=False)

def generare_list_with_structure(name,listfile='',structurefile=''):
    if not os.path.exists(listfile):
        list=pd.DataFrame()
        genes=f'working/{name}/{name}.fasta'
        records = SeqIO.parse(genes, 'fasta')
        for record in records:
            list = pd.concat([list, pd.DataFrame(
                {'Entry': [record.id.split('|')[1]], 'Entry Name': [record.id.split('|')[2]],
                 'Structure Data': [record.id.split('|')[1]]})])
        list.to_excel(f'working/{name}/{name}.xlsx', index=False)
    else:
        list=pd.read_excel(listfile)
        listdict=list.set_index('Structure Data')['Entry'].to_dict()
        for files in os.listdir(structurefile):
            try:
               os.rename(os.path.join(structurefile, files), os.path.join(structurefile, 'AF-'+listdict[files.split('.pdb')[0]]+'-F1-model_v4.pdb'))
               shutil.move( os.path.join(structurefile, 'AF-'+listdict[files]+'-F1-model_v4.pdb'),os.path.join(f'struct_data/taryeast/{name}', 'AF-'+listdict[files]+'-F1-model_v4.pdb'))
            except:
                pass
        list.to_excel(f'working/{name}/{name}.xlsx', index=False)

