import os
import pandas as pd
import shutil
import numpy as np
from Bio import SeqIO
def eggnog(name,refname,i=1):
    path=os.getcwd()
    tary = pd.read_excel(f'ziyuan/{name}.xlsx')
    homo = pd.read_excel(f'juzhen/juzhen_homolog{name}.xlsx')
    homo2 = []
    left = []
    for i in range(len(homo.index)):
        homo2.append(homo.iat[i, 2])

    for i in range(len(tary.index)):
        if tary.iat[i, 0] not in homo2:
            left.append(tary.iat[i, 0])
    left = pd.DataFrame(left)
    if i == 0:
        os.system('eggnog_mapper/download_eggnog_data.py')
        os.system('eggnog_mapper/create_dbs.py -m diamond --dbname yeast --taxa Saccharomycetes')
    if refname=='yeast':
        os.system(
        f'{path}/eggnog_mapper/emapper.py -m diamond -i {path}/ziyuan/{name}.fasta -o test{name} --tax_scope Saccharomycetes --cpu 0 --dmnd_db {path}/eggnog_mapper/data/yeast.dmnd --override')
    if refname=='ecoli':
        os.system(f'{path}/eggnog_mapper/emapper.py -m diamond -i {path}/ziyuan/{name}.fasta -o test{name} --tax_scope Bacteria --cpu 0 --dmnd_db {path}/eggnog_mapper/data/bacteria.dmnd --override')
    if refname=='strco':
        os.system(f'{path}/eggnog_mapper/emapper.py -m diamond -i {path}/ziyuan/{name}.fasta -o test{name} --tax_scope Actinobacteria --cpu 0 --dmnd_db {path}/eggnog_mapper/data/bacteria.dmnd --override')
    if refname=='human':
        os.system(f'{path}/eggnog_mapper/emapper.py -m diamond -i {path}/ziyuan/{name}.fasta -o test{name} --tax_scope Mammalia --cpu 0 --dmnd_db {path}/eggnog_mapper/data/mammalia.dmnd --override')
    with open(f'test{name}.emapper.annotations') as f:
        content = f.read()
        contents = content.split('\n')
    genes = pd.DataFrame()
    for i in contents:
        j = i.split('\t')
        j = np.transpose(pd.DataFrame((np.array(j))))
        genes = pd.concat([genes, j])
    genes.index = range(len(genes.index))
    for i in range(len(genes.index)):
        try:
            genes.iat[i, 0] = genes.iat[i, 0].split('|')[1]
        except:
            continue
    genes1 = []
    for i in range(len(genes.index)):
        genes1.append(genes.iat[i, 0])
    left.insert(loc=1, column='EC', value='')
    left.insert(loc=2, column='Rec', value='')
    for i in range(len(left.index)):
        try:
            left.iat[i, 1] = genes.iat[genes1.index(left.iat[i, 0]), 10]
            left.iat[i, 2] = genes.iat[genes1.index(left.iat[i, 0]), 14]
        except:
            left.iat[i, 1], left.iat[i, 2] = '-', '-'
    eggnogannop = pd.DataFrame()
    for i in range(len(left.index)):
        if left.iat[i, 2] != '-':
            eggnogannop = pd.concat([eggnogannop, left.iloc[[i], [0, 1, 2]]])
    eggnogannop.index = range(len(eggnogannop.index))
    eggnogannop.to_excel(f'juzhen/eggec2{name}.xlsx')
