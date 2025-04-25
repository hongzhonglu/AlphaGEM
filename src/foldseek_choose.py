import os
import pandas as pd
from plddt_find import get_plddt_from_pdb as gpl
from tqdm import tqdm
from multiprocessing import Pool

def process_file(path, filename):
    try:
        return filename.split('-')[1], gpl(os.path.join(path, filename))
    except:
        return filename.split('-')[1], None

def generate_dict(path, refname, pathwd):
    tardict = {}
    refdict = {}

    with Pool() as pool:
        results = list(tqdm(pool.starmap(process_file, [(path, i) for i in os.listdir(path)]), desc='generate dict1'))
        for result in results:
            tardict[result[0]] = result[1]

    with Pool() as pool:
        results = list(tqdm(pool.starmap(process_file, [(f'{pathwd}/struct_data/{refname}', i) for i in os.listdir(f'{pathwd}/struct_data/{refname}')]), desc='generate dict2'))
        for result in results:
            refdict[result[0]] = result[1]

    return tardict, refdict


def read_gx(name):
    global gx
    gx = pd.read_excel(f'./working/{name}/matrix_foldseek_{name}.xlsx')
gx2=pd.DataFrame()
def plddt_datacreate(name,refname,path=''):
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    tardict, refdict = generate_dict(path, refname, pathwd)
    return tardict,refdict



def tdblast(cmd1, cmd2, i, name, refname, path='', tardict=None, refdict=None):
    if tardict is None:
        tardict = {}
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    global gx2
    global gx
    # if gx.iat[i,5]<=0.8 or gx.iat[i,6]<=0.8:
    #     return 0
    d=tardict[cmd2]
    c=refdict[cmd1]
    a1=gx.iat[i,1]
    b1=gx.iat[i,2]
    a2=gx.iat[i,3]
    b2=gx.iat[i,4]
    a3=gx.iat[i,5]
    b3=gx.iat[i,6]
    gx2=pd.concat([gx2,pd.DataFrame({
        0:[a1],
        1:[b1],
        2:[a2],
        3:[b2],
        4:[a3],
        5:[b3],
        6:[c],
        7:[d]
    })])
def foldseek_choose(name,refname,path=''):
    read_gx(name)
    tardict,refdict=plddt_datacreate(name,refname)
    for i in tqdm(range(len(gx.index)),'find pLDDT'):
        tdblast(gx.iat[i, 1], gx.iat[i, 2], i, name,refname, path,tardict,refdict)
    print(gx2)
    gx2.to_excel(f'./working/{name}/matrix_foldseek_filtered1_{name}.xlsx')
