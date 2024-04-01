import pandas as pd
import numpy as np
import model_transpoter_found
gx = pd.read_excel('juzhen/juzhen4+ee.xlsx')
gx2 = pd.DataFrame()
gx3=pd.DataFrame()
count=0
for i in range(len(gx.index)):
    if model_transpoter_found.transpoter(gx.iat[i, 1]) == 1:
        if gx.iat[i, 5] >= 0.9 and gx.iat[i, 6] > 0.9 and gx.iat[i, 7] > 70 and gx.iat[i, 8] > 70 and gx.iat[
            i, 4] >= 0.5:
            gx2 = pd.concat([gx2, pd.DataFrame({
                0: [gx.iat[i, 1]],
                1: [gx.iat[i, 2]],
                2:[1],
                3:[gx.iat[i,4]]
            })])
    elif gx.iat[i, 5] >= 0.9 and gx.iat[i, 6] > 0.9 and gx.iat[i, 7] > 70 and gx.iat[i, 8] > 70 and gx.iat[
        i, 4] >= 0.5:
        gx2 = pd.concat([gx2, pd.DataFrame({
            0: [gx.iat[i, 1]],
            1: [gx.iat[i, 2]],
            2: [0],
            3:[gx.iat[i,4]]
        })])
gx2.columns=['names','tnames','ornot','tms']
yea = pd.read_excel('ziyuan/yeast.xlsx')
yea1 = []
yea2 = []
for i in range(len(yea.index)):
    yea1.append(yea.iat[i, 0])
    yea2.append(yea.iat[i, 2])
for i in range(len(gx2.index)):
    index = yea1.index(gx2.iat[i, 0])
    gx2.iat[i, 0] = yea2[index]
gx3=gx2.groupby('names').agg({'tms':['count','mean'],'ornot':['mean']})
gx3.to_excel('juzhen/whytarnsportor.xlsx')