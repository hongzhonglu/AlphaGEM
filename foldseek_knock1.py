
import pandas as pd
import os
import shutil
import time
df_foldseek_find = pd.DataFrame()
def foldseekfind(path_taryeast_structure='',name='',refname=''):
    global df_foldseek_find
    pathwd=os.getcwd()
    if path_taryeast_structure=='':
        path_taryeast_structure=f'{pathwd}/struct_data/taryeast/{name}'
    start_time = time.time()
    yea = pd.read_excel(f'data_available/{refname}.xlsx')
    try:
        os.mkdir(f'{pathwd}/struct_data/taryeast/{name}db')
    except:
        print('has created')
    os.system(
        f'{pathwd}/tools/foldseek/bin/foldseek createdb {path_taryeast_structure}/ {pathwd}/struct_data/taryeast/{name}db/{name}db')

    def bidui():
        global df_foldseek_find
        os.system(
            f'{pathwd}/tools/foldseek/bin/foldseek easy-search {pathwd}/struct_data/{refname} {pathwd}/struct_data/taryeast/{name}db/{name}db {pathwd}/data/aln{name}.csv tmpFolder --format-output "query,target,alntmscore,prob,qcov,tcov" --alignment-type 2 -e 0.00001')
        aln = pd.read_csv(f'{pathwd}/data/aln{name}.csv',sep='\t',names=['ref','tar','tms','pro','rcov','tcov'])
        aln2 = aln[aln['pro'] >= 1]
        aln2['ref']=aln2['ref'].apply(lambda x:x[3:x.rfind('-F1-')])
        aln2['tar']=aln2['tar'].apply(lambda x: x[3:x.rfind('-F1-')])
        aln2.loc[:,'pro']=aln2['tms']
        aln2.loc[:,'tms'] =1
        df_foldseek_find=aln2
        return 0

    def duiying(cmd):
        for i in range(len(yea.index)):
            if cmd == yea.iat[i, 2]:
                return yea.iat[i, 0]
        return 0

    shutil.copy(f'{pathwd}/data/alno/aln.csv', f'{pathwd}/data/aln{name}.csv')
    bidui()
    df_foldseek_find.to_excel(f'working/{name}/matrix_foldseek_{name}.xlsx')
    stop_time = time.time()
    print('the time that cost is' + str(stop_time - start_time))