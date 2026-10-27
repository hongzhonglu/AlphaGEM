import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import shutil
import subprocess
import sys
import runpy


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
    subprocess.run(
        [sys.executable, os.path.abspath(__file__), "clean-inference", f"{name}_homoleft"],
        cwd=clean_path,
        check=True,
    )
    shutil.move(f'{clean_path}/results/inputs/{name}_homoleft_maxsep.csv', f'working/{name}/{name}_homoleft_maxsep.csv')


def _run_clean_inference(fasta_data):
    """Keep CLEAN's ESM child in the same interpreter and propagate failures."""
    sys.path.insert(0, os.path.abspath("src"))
    from CLEAN import utils

    def extract_esm(fasta_name):
        subprocess.run(
            [sys.executable, os.path.abspath(__file__), "esm-extraction",
             f"data/{fasta_name}.fasta"],
            check=True,
        )

    utils.retrive_esm1b_embedding = extract_esm
    sys.argv = ["CLEAN_infer_fasta.py", "--fasta_data", fasta_data]
    runpy.run_path("CLEAN_infer_fasta.py", run_name="__main__")


def _run_esm_extraction(fasta_file):
    """Allow only the Namespace stored in ESM1b's legacy checkpoint."""
    import argparse
    from contextlib import nullcontext
    import torch

    sys.argv = ["esm/scripts/extract.py", "esm1b_t33_650M_UR50S",
                fasta_file, "data/esm_data", "--include", "mean"]
    # PyTorch >= 2.6 loads weights only by default; older releases lack this API.
    context = (torch.serialization.safe_globals([argparse.Namespace])
               if hasattr(torch.serialization, "safe_globals") else nullcontext())
    with context:
        runpy.run_path("esm/scripts/extract.py", run_name="__main__")
    

def clean_result(name):
    get_fasta(name)
    clean(name)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Internal CLEAN subprocess launcher")
    parser.add_argument("mode", choices=("clean-inference", "esm-extraction"))
    parser.add_argument("input")
    args = parser.parse_args()
    if args.mode == "clean-inference":
        _run_clean_inference(args.input)
    else:
        _run_esm_extraction(args.input)
