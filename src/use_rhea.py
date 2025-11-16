import os
import pandas as pd
import shutil
import numpy as np
from Bio import SeqIO

from plmsearch import rhea_search
from ss_predict_choose import ss_predict_choose_rhea

def process_rhea(name,refname,threshold,left):
    header=["ID","Query"]
    rhea_homo=pd.read_excel(f"./working/{name}/{name}_rhea_homolog{str(threshold)}.xlsx",index_col=0,header=None,names=header)
    rhea_homo.reset_index(drop=True,inplace=True)
    if rhea_homo.empty:
        print("No rhea homologs found")
        rhea_homo=pd.DataFrame({"ID":[],"Query":[]})
    rhea_homo = rhea_homo.drop(0)
    rhea_homo = rhea_homo[~rhea_homo["Query"].isin(left)]
    gene_rhea=pd.read_csv(f"./rhea/rhea2uniprot_sprot.tsv",sep="\t")
    gene2rheaid=gene_rhea.set_index("ID")["RHEA_ID"].to_dict()
    rhea_homo["Rhea_id"]=rhea_homo["ID"].map(gene2rheaid)
    rhea_homo.to_excel(f"./working/{name}/{name}_rheaid_plmsearch.xlsx")
def rhea(name,refname,threshold=0.9):
    path=os.getcwd()
    tary = []
    homo = pd.read_excel(f'./working/{name}/matrix_homolog{name}.xlsx')
    homo2 = []
    left = []
    
    for record in SeqIO.parse(f"./working/{name}/{name}.fasta","fasta"):
        tary.append(record.id)

    for i in range(len(homo.index)):
        homo2.append(homo.iat[i, 2])

    for i in tary:
        if i not in homo2:
            left.append(i)
    left = pd.DataFrame(left)

    rhea_search(name,refname,threshold)
    ss_predict_choose_rhea(name,threshold)
    process_rhea(name,refname,threshold,left)

