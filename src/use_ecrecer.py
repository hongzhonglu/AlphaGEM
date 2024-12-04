import os
import numpy as np
import pandas as pd
from Bio import SeqIO
def get_fasta(name):
    fasta_file=open(f'ziyuan/{name}.fasta')
    fasta=SeqIO.parse(fasta_file,'fasta')
    path = os.getcwd()
    tary = pd.read_excel(f'ziyuan/{name}.xlsx')
    homo = pd.read_excel('juzhen/juzhen_homolog.xlsx')
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
    fasta_out=open(f'ziyuan/{name}_homoleft.fasta','w')
    SeqIO.write(records_left,fasta_out,'fasta')
    fasta_out.close()

def ecrecer(name):
    os.system()

def ecrecer_result(name):
    get_fasta(name)
    ecrecer(name)

