import pandas as pd
import re
import os
from plddt_find import get_plddt_from_pdb as gpl
from functools import partial
import multiprocessing
from tqdm import tqdm
import subprocess

def parse_usalign_metrics(output):
    tm_scores = re.findall(r"TM-score=\s*([\d\.]+)", output)
    aligned_match = re.search(r"Aligned length=\s*(\d+)", output)
    length_matches = re.findall(
        r"TM-score=\s*[\d\.]+\s+\(normalized by length of Structure_[12]: L=(\d+),",
        output
    )

    if len(tm_scores) < 2:
        return '', '', ''

    tm_score_1 = float(tm_scores[0])
    tm_score_2 = float(tm_scores[1])

    if aligned_match and len(length_matches) >= 2:
        aligned_length = int(aligned_match.group(1))
        len_1 = int(length_matches[0])
        len_2 = int(length_matches[1])
        avg_coverage = round(((aligned_length / len_1) + (aligned_length / len_2)) / 2, 4)
    else:
        avg_coverage = ''


    return tm_score_1, tm_score_2, avg_coverage


def worker(i, gx, yeadict, name, refname, path):
    if gx.iat[i, 1]=='b448':
        print(duiying1(gx.iat[i, 1], yeadict))
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
    avg_cov=''
    pdd=pd.DataFrame({1:[0],2:[0],3:[0],4:[0],5:[0],6:[0],7:[0]})
    if cmd1==0 or cmd2==0:
        return pdd

    usalign_exec = os.path.join(pathwd, "tools/USalign/USalign")
    input_file_1 = os.path.join(path, f"AF-{cmd2}-F1-model_v4.pdb")
    input_file_2 = os.path.join(pathwd, f"struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb")

    cmd = [usalign_exec, input_file_1, input_file_2]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = res.stdout
        d, e, avg_cov = parse_usalign_metrics(output)

    except subprocess.CalledProcessError as e:
        print(f"Command execution failed: {e.stderr}")
        d = e = avg_cov = ''
    except Exception as ex:
        print(f"Parsing error occurred: {ex}")
        d = e = avg_cov = ''

    # print(f"Score 1: {d}, Score 2: {e}")
    # print('USalign done', cmd1, cmd2)
    if d!='':
            f=gpl(
            f'{path}/AF-{cmd2}-F1-model_v4.pdb')
            g=gpl(
            f'{pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb')
    else:
            f=''
            g=''
    pdd=pd.concat([
        pd.DataFrame([gx.iat[i,1]]),
        pd.DataFrame([gx.iat[i,2]]),
        pd.DataFrame([d]),
        pd.DataFrame([e]),
        pd.DataFrame([g]),
        pd.DataFrame([f]),
        pd.DataFrame([avg_cov]),
    ],axis=1)
    pdd.columns=[1,2,3,4,5,6,7]
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
