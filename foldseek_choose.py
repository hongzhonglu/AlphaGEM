import os
import pandas as pd
from plddt_find import get_plddt_from_pdb as gpl
def read_gx(name):
    global gx
    gx = pd.read_excel(f'juzhen/juzhen4{name}.xlsx')
gx2=pd.DataFrame()

def tdblast(cmd1,cmd2,i,name,refname,path=''):
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    global gx2
    global gx
    if gx.iat[i,5]<=0.8 or gx.iat[i,6]<=0.8:
        return 0
    d=gpl(f'{path}/AF-{cmd2}-F1-model_v4.pdb')
    c=gpl(
    f'{pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb')
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
    for i in range(len(gx.index)):
        tdblast(gx.iat[i, 1], gx.iat[i, 2], i, name,refname, path)
    print(gx2)
    gx2.to_excel(f'juzhen/juzhen4+ee{name}.xlsx')
