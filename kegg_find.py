import pandas as pd
from bioservices import KEGG
def find_tarname(name,strr):
    for i in strr:
        if name in i:
            return i[12:]
    return 'none'

def reaction_get():
    K = KEGG()
    reactions=[]
    reaction=pd.read_excel('juzhen/eggec2.xlsx')
    for i in range(len(reaction.index)):
         reactions+=reaction.iat[i,3].split(',')
    reactions=list(set(reactions))
    reactions2=pd.DataFrame(reactions,columns=['REACTION'])
    reactions2.insert(loc=1,column='NAME',value='')
    reactions2.insert(loc=2,column='DEFINITION',value='')
    reactions2.insert(loc=3,column='EQUATION',value='')
    n2=''
    n1=1
    while n1!=0:
        n1=0
        try:
            for i, r in reactions2.iterrows():
                r = r['REACTION']
                if reactions2.iat[i, 1] != '':
                    continue
                response = K.get(f'{r}')
                try:
                    rec2 = response
                    rec2 = rec2.split('\n')
                    reactions2.iat[i, 1] = find_tarname('NAME        ', rec2)
                    reactions2.iat[i, 2] = find_tarname('DEFINITION  ', rec2)
                    reactions2.iat[i, 3] = find_tarname('EQUATION    ', rec2)
                except:
                    n1 += 1
                    print('404 not find reaction {}'.format(r))
        except:
            continue
        if n2==n1:
            break
        n2=n1
    reactions2.to_excel('juzhen/reactions_kegg.xlsx')
def genes_get(name):
    reactions = []
    reaction = pd.read_excel('juzhen/eggec2.xlsx')
    for i in range(len(reaction.index)):
        reactions += reaction.iat[i, 3].split(',')
    reactions = list(set(reactions))
    gene_reaction_rule=pd.DataFrame()
    query = pd.read_excel(f'ziyuan/{name}.xlsx')
    for index,row in reaction.iterrows():
        for i in row['Rec'].split(','):
            gene_reaction_rule=pd.concat([gene_reaction_rule,pd.DataFrame({
                'genes':query.iat[list(query['Entry Name']).index(row[0]),0],
                'reaction':i
            },index=[0])])
    gene_reaction_rule2=gene_reaction_rule.groupby('reaction')['genes'].apply(lambda x:' or '.join(x)).reset_index()
    gene_reaction_rule2.to_excel('juzhen/gpr.xlsx')






