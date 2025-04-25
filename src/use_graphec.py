import os
import sys
import subprocess
import shutil
from Bio import SeqIO
def use_graphec(name):
    currentpath=os.getcwd()
    structures=[struc.split('-')[1] for struc in os.listdir(f'./struct_data/taryeast/{name}')]
    genes=SeqIO.parse(f'working/{name}/{name}_homoleft.fasta','fasta')

    os.chdir('./tools/GraphEC/')
    os.system(f"python main.py --task EC_number --fasta ../working/{name}/{name}_homoleft.fasta --name {name}")
    shutil.move('./EC_number/results/example_top5.txt',f'../working/{name}/{name}graphecresult.csv')
    os.chdir(currentpath)