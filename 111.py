import re
import re

import pandas as pd


def distribute_or_over_and(expr):
    # 正则表达式匹配'or'和'and'以及括号
    pattern = r'\(([^()]+)\)'
    or_pattern = r'(.+?) or (.+?)'

    # 递归处理括号内的表达式
    while '(' in expr and ')' in expr:
        expr = re.sub(pattern, lambda m: distribute_or_over_and(m.group(1)), expr)

    # 处理'or'运算符，应用分配律
    while re.search(or_pattern, expr):
        expr = re.sub(or_pattern, lambda m: f'({m.group(1)}) or ({m.group(2)})', expr)

    # 处理'and'运算符，确保'or'提升到外层
    expr = expr.replace(' and ', ') and (')

    return expr


# 示例逻辑表达式
logic_expr = '(A and B) or (C and (D or E)) or (F and G)'

# 应用分配律
new_logic_expr = distribute_or_over_and(logic_expr)
print(new_logic_expr)


df=pd.DataFrame()
for genes in model.genes:
    try:
        df=pd.concat([df,pd.DataFrame({'genename':genes.id,'uniprot':genes.annotation.get('uniprot')},index=[1])])
    except:
        df = pd.concat(
            [df, pd.DataFrame({'genename': genes.id, 'uniprot': genes.annotation.get('uniprot')[0]}, index=[1])])