import pandas as pd
import numpy as np
def US_align_choose(name):
    gx3 = pd.DataFrame()
    gx2 = pd.read_excel(f'./working/{name}/matrix_USalign_filtered{name}.xlsx')
    gx2.fillna(value=1, inplace=True)
    a = 0
    for i in range(len(gx2.index)):
        if gx2.iat[i, 3] + gx2.iat[i, 4] >= 1.4 or gx2.iat[i, 5] + gx2.iat[i, 6] <= 1.6:
            gx3 = pd.concat([gx3, pd.DataFrame({0: [gx2.iat[i, 1]],
                                                1: [gx2.iat[i, 2]],
                                                2: [gx2.iat[i, 3]],
                                                3: [gx2.iat[i, 4]],
                                                4: [gx2.iat[i, 5]],
                                                5: [gx2.iat[i, 6]]
                                                })])
            a += 1
    gx3.to_excel(f'working/{name}/matrix_homo_part1{name}.xlsx')