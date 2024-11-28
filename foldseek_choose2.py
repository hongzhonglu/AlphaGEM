import pandas as pd
import numpy as np
import model_transpoter_found
import align2
def foldseek_choose(refname,refmodel,name,TMscoretrans,coveragetrans,TMscore,coverage,pLDDt):
    gx = pd.read_excel(f'working/{name}/matrix_foldseek_filtered1_{name}.xlsx')
    model_transpoter_found.datapre(refmodel,refname)
    align2.predate(name,refname)
    gx2 = pd.DataFrame()
    for i in range(len(gx.index)):
        if model_transpoter_found.transpoter(gx.iat[i, 1]) == 1:
            if gx.iat[i, 5] > TMscoretrans and gx.iat[i, 6] > TMscoretrans and gx.iat[i, 7] > pLDDt and gx.iat[i, 8] > pLDDt and gx.iat[i, 4] >= coveragetrans:
                gx2 = pd.concat([gx2, pd.DataFrame({
                    'gene1': [gx.iat[i, 1]],
                    'gene2': [gx.iat[i, 2]],
                    'coverage': [(gx.iat[i, 5]+gx.iat[i, 6])/2],
                    'tmscore': [gx.iat[i, 4]],
                    'sumscore': [(gx.iat[i, 5]+gx.iat[i, 6])/2+gx.iat[i, 4]],
                    'transporter':[True]
                })])
        elif gx.iat[i, 5] > TMscore and gx.iat[i, 6] > TMscore and gx.iat[i, 7] > pLDDt and gx.iat[i, 8] > pLDDt and gx.iat[
            i, 4] >= coverage:
            gx2 = pd.concat([gx2, pd.DataFrame({
                'gene1': [gx.iat[i, 1]],
                'gene2': [gx.iat[i, 2]],
                'coverage': [(gx.iat[i, 5] + gx.iat[i, 6])/2],
                'tmscore': [gx.iat[i, 4]],
                'sumscore': [(gx.iat[i, 5] + gx.iat[i, 6]) / 2 + gx.iat[i, 4]],
                'transporter': [False]
            })])
    yea = pd.read_excel(f'data_available/{refname}.xlsx')
    yea1 = []
    yea2 = []
    for i in range(len(yea.index)):
        yea1.append(yea.iat[i, 0])
        yea2.append(yea.iat[i, 2])
    for i in range(len(gx2.index)):
        index = yea1.index(gx2.iat[i, 0])
        gx2.iat[i, 0] = yea2[index]
    gx2.to_excel(f'./working/{name}/matrix_foldseek_filtered2_{name}.xlsx')
def foldseek_choose2(refname,refmodel,name,TMscoretrans=0.7,coveragetrans=0.8,TMscore=0.7,coverage=0.8,pLDDt=70):
    foldseek_choose(refname,refmodel,name,coveragetrans,TMscoretrans,coverage,TMscore,pLDDt)