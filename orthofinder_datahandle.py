import numpy as np
import pandas as pd
import os
import shutil
def datahandel(name='',refname=''):
    orth = pd.DataFrame()
    path=os.getcwd()
    fasta=os.listdir(f'{path}/orth/data')
    for i in fasta[0:-1]:
        try:
          os.remove(os.path.join(path+'/orth/data',i))
        except:
            continue

    shutil.copy(f'ziyuan/{refname}.fasta', f'orth/data/z-{refname}.fasta')
    shutil.copy(f'ziyuan/{name}.fasta', f'orth/data/a-{name}.fasta')
    os.system(
        f'{path}/tools/OrthoFinder/orthofinder -f {path}/orth/data/')
    pathresult=f'{path}/orth/data/OrthoFinder/'+os.listdir(f'{path}/orth/data/OrthoFinder/')[0]+'/Orthogroups/Orthogroups.tsv'
    with open(
           pathresult) as orth1:  # key is the result of orthofinder
        orth2 = orth1.read()
        orth3 = orth2.split('\n')
        for i in orth3[1:]:
            orth = pd.concat([orth, np.transpose(pd.DataFrame(i.split('\t')))])
    tyea = pd.DataFrame()  # kid is the yeast while the klg and kid2 is the yil
    yea = pd.DataFrame()
    orth = orth.fillna('none')
    orth.astype(str)
    for i in range(len(orth.index)):
        tyea2 = pd.DataFrame(np.zeros([1, 200]))
        num = 0
        num2 = 0
        l = 0
        bb = 0
        for k in range(len(orth.iat[i, 1])):
            if orth.iat[i, 1][num] == '|':
                bb += 1
                if bb % 2 == 1:
                    num2 = num + 1
                if bb % 2 == 0:
                    tyea2.iat[0, l] = orth.iat[i, 1][num2:num]
                    l += 1
            num += 1
        tyea = pd.concat([tyea, tyea2])
    yea2 = pd.DataFrame()
    for i in range(len(orth.index)):
        yea3 = pd.DataFrame(np.zeros([1, 200]))
        num = 0
        num2 = 0
        l = 0
        bb = 0
        for k in range(len(orth.iat[i, 2])):
            if orth.iat[i, 2][num] == '|':
                bb += 1
                if bb % 2 == 1:
                    num2 = num + 1
                if bb % 2 == 0:
                    yea3.iat[0, l] = orth.iat[i, 2][num2:num]
                    l += 1
            num += 1
        yea2 = pd.concat([yea2, yea3])
    ref = pd.read_excel(
        f'{path}/ziyuan/{refname}.xlsx')
    yea = yea2
    juzhen1 = pd.DataFrame()
    for i in range(len(tyea.index)):
        for j in range(200):
            if tyea.iat[i, j] == 0.0:
                break
            for m in range(200):
                if yea2.iat[i, m] == 0.0:
                    break
                jp = pd.concat([pd.DataFrame([yea2.iat[i, m]]), pd.DataFrame([tyea.iat[i, j]])], axis=1)
                juzhen1 = pd.concat([juzhen1, jp])
    ref2=[]
    for i in range(len(ref.index)):
        ref2.append(ref.iat[i,0])
    for i in range(len(juzhen1.index)):
            try:
                juzhen1.iat[i,0]=ref.iat[ref2.index(juzhen1.iat[i,0]),2]
            except:
                continue
    print(yea)
    juzhen1.to_excel('juzhen/juzhen1.xlsx')
