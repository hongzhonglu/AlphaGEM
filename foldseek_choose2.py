import pandas as pd
import numpy as np
import model_transpoter_found
import align2
def foldseek_choose(refname):
    gx = pd.read_excel('juzhen/juzhen4+ee.xlsx')
    gx2 = pd.DataFrame()
    for i in range(len(gx.index)):
        if model_transpoter_found.transpoter(gx.iat[i, 1]) == 1:
            if gx.iat[i, 5] >= 0.9 and gx.iat[i, 6] > 0.9 and gx.iat[i, 7] > 70 and gx.iat[i, 8] > 70 and gx.iat[i, 4] >= 0.85 and align2.align(gx.iat[i,1],gx.iat[i,2])>=0.5:
                gx2 = pd.concat([gx2, pd.DataFrame({
                    0: [gx.iat[i, 1]],
                    1: [gx.iat[i, 2]]
                })])
        elif gx.iat[i, 5] >= 0.9 and gx.iat[i, 6] > 0.9 and gx.iat[i, 7] > 70 and gx.iat[i, 8] > 70 and gx.iat[
            i, 4] >= 0.5:
            gx2 = pd.concat([gx2, pd.DataFrame({
                0: [gx.iat[i, 1]],
                1: [gx.iat[i, 2]]
            })])
    yea = pd.read_excel(f'ziyuan/{refname}.xlsx')
    yea1 = []
    yea2 = []
    for i in range(len(yea.index)):
        yea1.append(yea.iat[i, 0])
        yea2.append(yea.iat[i, 2])
    for i in range(len(gx2.index)):
        index = yea1.index(gx2.iat[i, 0])
        gx2.iat[i, 0] = yea2[index]
    gx2.to_excel('juzhen/juzhen_model3.xlsx')
def foldseek_choose2(refname):
    foldseek_choose(refname)