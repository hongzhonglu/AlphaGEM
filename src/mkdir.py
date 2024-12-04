import os
import shutil
import pandas
import pandas as pd


def makeworkdir(name,fasta):
    directory_path=f"./working/{name}"
    try:
        os.mkdir(directory_path)
        print("Directory created")
    except FileExistsError:
        print("Dir exists")
    except OSError as e:
        print(f"Error:{e}")

    try:
        shutil.move(fasta,directory_path)
    except:
        print("file exists")


def mvfile(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction):
    df=pd.read_excel(f"./working/{name}/{name}_ss{threshold}_cov{str(cov_threshold)}_id{str(id_threshold)}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.xlsx")
    df=df[["Gene1","Gene2"]]
    df.columns=[0,1]
    df.to_excel(f"./working/{name}/matrix_homolog{name}.xlsx")