import pandas as pd

def homo(name):
    gx1 = pd.read_excel(f'working/{name}/matrix_homo_part1{name}.xlsx')
    gx4 = pd.read_excel(f'working/{name}/matrix_homo_part2{name}.xlsx')
    
    # 使用 set 追踪已添加的 (gene, homolog) 配对
    pairs = set()
    rows = []
    
    # 从 gx1 提取同源对
    for i in range(len(gx1.index)):
        gene = gx1.iat[i, 1]
        homolog = gx1.iat[i, 2]
        pair = (gene, homolog)
        if pair not in pairs:
            pairs.add(pair)
            rows.append({0: gene, 1: homolog})
    
    # 从 gx4 添加不重复的同源对
    for i in range(len(gx4.index)):
        gene = gx4.iat[i, 1]
        homolog = gx4.iat[i, 2]
        pair = (gene, homolog)
        if pair not in pairs:
            pairs.add(pair)
            rows.append({0: gene, 1: homolog})
    
    gx3 = pd.DataFrame(rows)
    gx3.to_excel(f'working/{name}/matrix_homolog{name}_preforprottrans.xlsx', index=False)