import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import shutil
import subprocess
import sys


def activate_clean_env():
    try:
        subprocess.run("source activate clean", shell=True, check=True)
        print("Switched to CLEAN environment.")
    except subprocess.CalledProcessError as e:
        print(f"Error switching to CLEAN environment: {e}")


def activate_alphagem_env():
    try:
        subprocess.run("source activate AlphaGEM", shell=True, check=True)
        print("Switched to AlphaGEM environment.")
    except subprocess.CalledProcessError as e:
        print(f"Error switching to AlphaGEM environment: {e}")
path=os.getcwd()
clean_path=os.getcwd()+'/tools/CLEAN/app'
def get_fasta(name):
    fasta_file=open(f'working/{name}/{name}.fasta')
    fasta=SeqIO.parse(fasta_file,'fasta')
    path = os.getcwd()
    tary = pd.read_excel(f'working/{name}/{name}.xlsx')
    homo = pd.read_excel(f'working/{name}/matrix_homolog{name}.xlsx')
    homo2 = []
    left = []
    records_left=[]
    for i in range(len(homo.index)):
        homo2.append(homo.iat[i, 2])
    for i in range(len(tary.index)):
        if tary.iat[i, 0] not in homo2:
            left.append(tary.iat[i, 0])
    for records in fasta:
        if records.id.split('|')[1] in left:
            records_left.append(records)
    fasta_out=open(f'{clean_path}/data/inputs/{name}_homoleft.fasta','w')
    SeqIO.write(records_left,fasta_out,'fasta')
    fasta_out.close()
    fasta_out = open(f'working/{name}/{name}_homoleft.fasta', 'w')
    SeqIO.write(records_left, fasta_out, 'fasta')
    fasta_out.close()

def clean(name):
        #activate_clean_env()
        os.chdir(f'{clean_path}')
        sys.path.append('src')
        os.system(f'python {clean_path}/CLEAN_infer_fasta.py --fasta_data {name}_homoleft')
        os.chdir(f'{path}')
        shutil.move(f'{clean_path}/results/inputs/{name}_homoleft_maxsep.csv', f'working/{name}/{name}_homoleft_maxsep.csv')
        #activate_alphagem_env()
    

def clean_result(name):
    get_fasta(name)
    clean(name)

