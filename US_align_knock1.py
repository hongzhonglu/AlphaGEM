import pandas as pd
import os
from plddt_find import get_plddt_from_pdb as gpl
def read_gx():
    global gx
    gx = pd.read_excel('juzhen/juzhen1.xlsx')

def tdblast(cmd1,cmd2,i,name,refname,path=''):
    pathwd=os.getcwd()
    if path=='':
        path=f'{pathwd}/struct_data/taryeast/{name}'
    d=1
    e=1
    pdd=pd.DataFrame({1:[0],2:[0],3:[0],4:[0],5:[0],6:[0]})
    if cmd1==0 or cmd2==0:
        return pdd
    res = os.popen(f"{pathwd}/tools/USalign/USalign {path}/AF-{cmd2}-F1-model_v4.pdb {pathwd}/struct_data/{refname}/AF-{cmd1}-F1-model_v4.pdb")
    for j in range(20):
        r = res.readline()
        print(r)
        try:
           if j==10:
               a=eval(r[23:26])
           if j==11:
               b=eval(r[23:26])
           if j==13:
               c=eval(r[16:19])
           if j == 14:
               d=r[10:16]
           if j == 15:
               e=r[10:16]
        except:
            d=''
            e=''
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
    for i in range(len(yea.index)):
        if a==yea.iat[i,2]:
            return yea.iat[i,0]
    return 0
def duiying2(a,tar):
    for i in range(len(tar.index)):
        if a==tar.iat[i,1]:
            return tar.iat[i,0]
    return 0

def US_align_find(name,path,refname,path2):
    yea = pd.read_excel(f'ziyuan/{refname}.xlsx')
    read_gx()
    gx2 = pd.DataFrame()
    tar = pd.read_excel(f'ziyuan/{name}.xlsx')
    for i in range(len(gx.index)):
        gx2 = pd.concat([gx2, tdblast(duiying1(gx.iat[i, 1],yea), gx.iat[i, 2], i,name,refname,path)])
    gx2.index=range(len(gx2.index))
    print(gx2)
    gx2.to_excel('juzhen/juzhen2.xlsx')