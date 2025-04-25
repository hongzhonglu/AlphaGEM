import pandas as pd
import os
from plddt_find import get_plddt_from_pdb as gpl
from functools import partial
import multiprocessing
from tqdm import tqdm
import subprocess


def worker(i, gx, yeadict, name, refname, path):
    return tdblast(duiying1(gx.iat[i, 1], yeadict), gx.iat[i, 2], i, name, refname, path)


def parallel_processing(gx, yeadict, name, refname, path):
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    worker_with_args = partial(worker, gx=gx, yeadict=yeadict, name=name, refname=refname, path=path)
    results = list(tqdm(pool.imap(worker_with_args, range(len(gx.index))), total=len(gx.index)))
    pool.close()
    pool.join()

    return results

def read_gx(name):
    global gx
    gx = pd.read_excel(f'working/{name}/matrix_orthofinder{name}.xlsx')

def tdblast(cmd1,cmd2,i,name,refname,path=''):
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    d=1
    e=1
    pdd=pd.DataFrame({1:[0],2:[0],3:[0],4:[0],5:[0],6:[0]})
    if cmd1==0 or cmd2==0:
        return pdd
    cmd = f"{pathwd}/tools/USalign/USalign {path}/AF-{cmd2}-F1-model_v4.pdb {pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb"
    res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = res.stdout.readlines()
    try:
        d = float(output[15][10:17])
        e = float(output[16][10:17])
    except Exception as ex:
        print(f"Error occurred: {ex}")
        d = e = ''
    print('done', cmd1, cmd2)
    if d!='':
            f=gpl(
            f'{path}/AF-{cmd2}-F1-model_v4.pdb')
            g=gpl(
            f'{pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb')
    else:
            f=''
            g=''
    pdd=pd.concat([pd.DataFrame([gx.iat[i,1]]),pd.DataFrame([gx.iat[i,2]]),pd.DataFrame([d]),pd.DataFrame([e]),pd.DataFrame([g]),pd.DataFrame([f])],axis=1)
    pdd.columns=[1,2,3,4,5,6]
    return pdd

def duiying1(a,yea):
    try:
        return yea[a]
    except:
        return a
def duiying2(a,tar):
    for i in range(len(tar.index)):
        if a==tar.iat[i,1]:
            return tar.iat[i,0]
    return 0

def US_align_find(name,path,refname,path2):
    yea = pd.read_excel(f'data_available/{refname}.xlsx')
    yeadict={yea.iat[i,2]:yea.iat[i,0] for i in range(len(yea.index))}
    read_gx(name)
    gx2 = pd.DataFrame()
    tar = pd.read_excel(f'working/{name}/{name}.xlsx')
    results=[]
    print('start usalign alignment')
    results = parallel_processing(gx, yeadict, name, refname, path)
    gx2=pd.concat(results)
    gx2.index=range(len(gx2.index))
    print(gx2)
    gx2.to_excel(f'working/{name}/matrix_USalign_filtered{name}.xlsx')