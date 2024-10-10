import os
from use_clean import activate_alphagem_env
import subprocess

def activate_deepectransformer_env():
    try:
        subprocess.run("source activate deepectransformer", shell=True, check=True)
        print("Switched to deepectransformer environment.")
    except subprocess.CalledProcessError as e:
        print(f"Error switching to CLEAN environment: {e}")


def use_deepectransformer(name):
    path = os.getcwd()
    os.chdir(f"{path}/DeepProZyme")
    #activate_deepectransformer_env()
    os.system(f'conda run -n deepectransformer python {path}/DeepProZyme/run_deepectransformer.py -i {path}/ziyuan/{name}_homoleft.fasta -o {path}/juzhen/{name}_deepec_result -g cpu -b 128 -cpu 2')
    os.chdir(f"{path}")
    #activate_alphagem_env()