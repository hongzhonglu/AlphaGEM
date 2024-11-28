import gc
import cobra
from rdflib import Graph, Namespace, URIRef
import pandas as pd
import tqdm
import concurrent.futures
from multiprocessing import Manager
import pickle
from multicpu import multi_cpu
def numzzp(x):
    return [x['metacyc'], x['kegg'], x['EC'], x['GO']]
rhea2databases = pd.read_excel('data_available/rhea2database.xlsx')
rhea2databases = rhea2databases.fillna('')
rhea2datadict = rhea2databases.set_index(['ID']).apply(numzzp, axis=1).to_dict()
rhea2direction = rhea2databases.set_index(['ID'])['direction'].to_dict()
kegg2namepd = pd.read_csv('data_available/reaction', sep='\t', header=None, names=['id', 'name'])
kegg2name = kegg2namepd.set_index('id')['name'].to_dict()


def add_reactions(resultss):
    reactions,x,reaction_dict,compoundsdict,genericcpdlist,chebicpdlist=resultss
    result_dict = {}
    reaction_add = cobra.Reaction()
    gc.collect()
    for reaction in reactions:
      try:
        contains_list = reaction_dict[reaction]
        for contain, contains in contains_list:
            row=compoundsdict[contain]
            metabolite = cobra.Metabolite()
            met = str(row.a)
            if met.startswith('GENERIC'):
                row = genericcpdlist[met]
                metabolite = cobra.Metabolite(
                    name=f"{row.n}",
                    charge=int(f"{row.c}"),
                    id=f"{met}",
                    formula=f"{row.f}",
                    compartment='c'
                )
                metabolite.annotation['chebi'] = str(row.chebi).split('_')[1]
            if met=='CHEBI:10545':
                continue
            if met.startswith('CHEBI'):
                row = chebicpdlist[met]
                metabolite = cobra.Metabolite(
                    name=f"{row.n}",
                    charge=int(f"{row.c}"),
                    id=f"{met}",
                    formula=f"{row.f}",
                    compartment='c'
                )
                metabolite.annotation['chebi'] = str(row.chebi).split('_')[1]
            if contains.split('/')[-1][8:] == 'N' or contains.split('/')[-1][8:] == 'n':
                if reaction.endswith('L'):
                    result_dict[metabolite] = -4
                else:
                    result_dict[metabolite] = 4
            elif contains.split('/')[-1][8:] == 'Nplus1':
                result_dict[metabolite] = -5 if x.endswith('L') else 5
            elif contains.split('/')[-1][8:] == 'Nminus1':
                result_dict[metabolite] = -3 if x.endswith('L') else 3
            elif contains.split('/')[-1][8:] == '2n':
                if reaction.endswith('L'):
                    result_dict[metabolite] = -8
                else:
                    result_dict[metabolite] = 8
            else:
                if reaction.endswith('L'):
                    result_dict[metabolite] = -int(contains.split('/')[-1][8:])
                else:
                    result_dict[metabolite] = int(contains.split('/')[-1][8:])
      except:
        pass

    reaction_add.lower_bound = 0
    reaction_add.upper_bound = 1000
    reaction_add.id = f'rhea_{x.split("/")[-1]}'
    reaction_add.annotation['rhea'] = x.split("/")[-1]
    reaction_add = add_anno(reaction_add, x)
    reaction_add.add_metabolites(result_dict)
    print(reaction_add.id)
    return reaction_add


# reaction的父子 关系，主要加入子反应。
# reaction以及metabolites的标准化。


def add_anno(reaction_add, x,):
    # 添加注释信息
    try:
        reaction_add.annotation['kegg'] = rhea2datadict[int(x.split("/")[-1])][1]
        reaction_add.name = kegg2name.get(reaction_add.annotation['kegg'], '').split(';')[0]
    except:
        pass

    keys = ['MetaCyc', 'GO', 'EC']
    for i, key in enumerate(keys):
        try:
            reaction_add.annotation[key] = rhea2datadict[int(x.split("/")[-1])][i]
        except:
            pass

    try:
        direction = rhea2direction[x.split("/")[-1]]
        if direction == 'RL':
            reaction_add.lower_bound, reaction_add.upper_bound = -1000, 0
        elif direction == 'BI':
            reaction_add.lower_bound, reaction_add.upper_bound = -1000, 1000
    except KeyError:
        pass

    return reaction_add


def main():
    manager = Manager()
    Reactions=[]
    g = Graph()
    # 从文件或URL中解析RDF数据
    g.parse('rhea.rdf', format='xml')
    rh = Namespace("http://rdf.rhea-db.org/")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    reaction_query = f"""
    SELECT DISTINCT ?x ?reaction ?general
    WHERE {{
      ?x rdfs:subClassOf rh:Reaction ;
         rdfs:subClassOf ?general ;
         rh:accession ?accession;
         rh:side ?reaction .
    }}
    """
    reactions = g.query(reaction_query)

    # 将结果存储在一个字典中
    reaction_dict = {str(row.reaction): [] for row in reactions}
    x_reaction_dict = {str(row.x): [] for row in reactions}
    for row in reactions:
        x_reaction_dict[str(row.x)].append(str(row.reaction))
    generaldict = {str(row.x).split('/')[-1]: [] for row in reactions}
    for row in reactions:
        if str(row.general).split('/')[-1] !='Reaction':
             generaldict[str(row.x).split('/')[-1]].append(str(row.general).split('/')[-1])
    for i in generaldict.keys():
        generaldict[i]=list(set(generaldict[i]))
    reactionsonsdict = {str(row.x).split('/')[-1]: [] for row in reactions}
    for i in generaldict.keys():
        for j in generaldict[i]:
            try:
                reactionsonsdict[j].append(i)
            except:
                print(j)
    nameskeys=list(reactionsonsdict.keys())
    for i in nameskeys:
        if len(reactionsonsdict[i]) ==0:
            del reactionsonsdict[i]
    fsave=open('data_available/general2specific.pkl', 'wb')
    pickle.dump(reactionsonsdict, fsave)
    fsave.close()
    fopen=open('data_available/general2specific.pkl', 'rb')
    a=pickle.load(fopen)
    fopen.close()
    # 第二步：基于第一步的结果，获取相关的 ?contain 和 ?contains
    contain_query = """
    SELECT ?reaction ?contain ?contains
    WHERE {
      VALUES ?reaction { """ + " ".join(f"<{r}>" for r in reaction_dict.keys()) + """ }
      ?reaction ?contains ?contain .
      ?contains rdfs:subPropertyOf rh:contains .
    }
    """
    contains = g.query(contain_query)

    # 将结果存储在字典中
    for row in contains:
        reaction_dict[str(row.reaction)].append((str(row.contain), str(row.contains)))

    # 第三步：基于第二步的结果，获取相关的化合物信息
    compound_query = """
    SELECT ?reaction ?contain ?contains ?a
    WHERE {
      VALUES ?contain { """ + " ".join(f"<{c[0]}>" for r in reaction_dict.values() for c in r) + """ }
      ?contain rh:compound ?cpd .
      ?cpd rh:accession ?a .
    }
    """

    # addreaction
    compounds = g.query(compound_query)
    compoundsdict = {str(row.contain): row for row in compounds}
    del compounds
    cpdquery = f"""
    SELECT DISTINCT ?cpd ?n ?c ?f ?chebi ?met
    WHERE {{
      ?cpd  rh:accession ?met ;
            rh:name ?n ;
            rh:reactivePart ?name .
      ?name rh:charge ?c ;
            rh:formula ?f ;
            rh:chebi ?chebi.}}
    """
    genericcpdlsit = g.query(cpdquery)
    genericcpdlist = {str(row.met): [] for row in genericcpdlsit}
    for row in genericcpdlsit:
        genericcpdlist[str(row.met)].append(row)
    del genericcpdlsit
    cpdquery = f"""
    SELECT DISTINCT ?cpd ?n ?c ?f ?met ?chebi
    WHERE {{
      ?cpd  rh:accession ?met ;
            rh:name ?n ;
            rh:charge ?c ;
            rh:formula ?f ;
            rh:chebi ?chebi.}}
    """
    chebicpdlsit = g.query(cpdquery)
    chebicpdlist = {str(row.met): row for row in chebicpdlsit}
    del chebicpdlsit
    cpdquery = f"""
        SELECT DISTINCT ?cpd ?n ?c ?f ?met ?chebi
        WHERE {{
          ?cpd  rh:accession ?met ;
                rh:name ?n ;
                rh:charge ?c ;
                rh:formula ?f ;
                rh:underlyingChebi ?chebi.}}
        """
    polymercpdlsit = g.query(cpdquery)
    polymercpdlist = {str(row.met): row for row in polymercpdlsit}
    del polymercpdlsit
    for x,reactions in tqdm.tqdm(x_reaction_dict.items()):
        result_dict = {}
        reaction_add = cobra.Reaction()
        for reaction in reactions:
                contains_list = reaction_dict[reaction]
                for contain, contains in contains_list:
                    row = compoundsdict[contain]
                    metabolite = cobra.Metabolite()
                    met = str(row.a)
                    if met.startswith('GENERIC'):
                      c=0
                      for row in genericcpdlist[met]:
                        c+=1
                        metabolite = cobra.Metabolite(
                            name=f"{row.n}",
                            charge=int(f"{row.c}"),
                            id=f"{met}"+'part'+str(c),
                            formula=f"{row.f}",
                            compartment='c'
                        )
                        metabolite.annotation['chebi'] = str(row.chebi).split('_')[1]
                        if contains.split('/')[-1][8:] == 'N' or contains.split('/')[-1][8:] == 'n':
                            if reaction.endswith('L'):
                                result_dict[metabolite] = -4
                            else:
                                result_dict[metabolite] = 4
                        elif contains.split('/')[-1][8:] == 'Nplus1':
                            result_dict[metabolite] = -5 if x.endswith('L') else 5
                        elif contains.split('/')[-1][8:] == 'Nminus1':
                            result_dict[metabolite] = -3 if x.endswith('L') else 3
                        elif contains.split('/')[-1][8:] == '2n':
                            if reaction.endswith('L'):
                                result_dict[metabolite] = -8
                            else:
                                result_dict[metabolite] = 8
                        else:
                            if reaction.endswith('L'):
                                result_dict[metabolite] = -int(contains.split('/')[-1][8:])
                            else:
                                result_dict[metabolite] = int(contains.split('/')[-1][8:])
                      continue
                    elif met == 'CHEBI:10545' or met == 'CHEBI:30212':
                        continue
                    elif met.startswith('POLYMER'):
                        row = polymercpdlist[met]
                        metabolite = cobra.Metabolite(
                            name=f"{row.n}",
                            charge=int(f"{row.c}".split(')(')[0][1:]) if '(' in f"{row.c}" else int(f"{row.c}"),
                            id=f"{met}",
                            formula=f"{row.f}".split('(')[0],
                            compartment='c'
                        )
                        metabolite.annotation['chebi'] = str(row.chebi).split('_')[1]
                    else:
                        row = chebicpdlist[met]
                        metabolite = cobra.Metabolite(
                            name=f"{row.n}",
                            charge=int(f"{row.c}"),
                            id=f"{met}",
                            formula=f"{row.f}",
                            compartment='c'
                        )
                        metabolite.annotation['chebi'] = str(row.chebi).split('_')[1]
                    if contains.split('/')[-1][8:] == 'N' or contains.split('/')[-1][8:] == 'n':
                        if reaction.endswith('L'):
                            result_dict[metabolite] = -4
                        else:
                            result_dict[metabolite] = 4
                    elif contains.split('/')[-1][8:] == 'Nplus1':
                        result_dict[metabolite] = -5 if x.endswith('L') else 5
                    elif contains.split('/')[-1][8:] == 'Nminus1':
                        result_dict[metabolite] = -3 if x.endswith('L') else 3
                    elif contains.split('/')[-1][8:] == '2n':
                        if reaction.endswith('L'):
                            result_dict[metabolite] = -8
                        else:
                            result_dict[metabolite] = 8
                    else:
                        if reaction.endswith('L'):
                            result_dict[metabolite] = -int(contains.split('/')[-1][8:])
                        else:
                            result_dict[metabolite] = int(contains.split('/')[-1][8:])


        reaction_add.lower_bound = 0
        reaction_add.upper_bound = 1000
        reaction_add.id = f'rhea_{x.split("/")[-1]}'
        reaction_add.annotation['rhea'] = x.split("/")[-1]
        reaction_add = add_anno(reaction_add, x)
        reaction_add.add_metabolites(result_dict)
        Reactions.append(reaction_add)
    #Reactions=multi_cpu(add_reactions,[[reactions,x,reaction_dict,compoundsdict,genericcpdlist,chebicpdlist] for x, reactions in tqdm.tqdm(x_reaction_dict.items())],8,2)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # 提交任务给线程池
        futures = [
            executor.submit(add_reactions, [reactions, x, reaction_dict, compoundsdict, genericcpdlist, chebicpdlist])
            for x, reactions in tqdm.tqdm(x_reaction_dict.items(), desc="Processing Reactions")]
        # 获取结果
        for future in concurrent.futures.as_completed(futures):
            Reactions.append(future.result())
    model = cobra.Model('universial')
    model.add_reactions(Reactions)
    for r in tqdm.tqdm(Reactions):
      try:
        model.add_reactions([r])
      except Exception as e:
        print(e)
    cobra.io.save_json_model(model, 'models/reaction_pool.json')

main()


# import cobra
# from rdflib import Graph, Namespace, URIRef
# import pandas as pd
# from tqdm import tqdm
#
#
#
#
# # 读取和处理数据
# rhea2databases = pd.read_excel('data_available/rhea2database.xlsx').fillna('')
# rhea2datadict = rhea2databases.set_index(['ID']).apply(numzzp, axis=1).to_dict()
# rhea2direction = rhea2databases.set_index(['ID'])['direction'].to_dict()
# kegg2namepd = pd.read_csv('data_available/reaction', sep='\t', header=None, names=['id', 'name'])
# kegg2name = kegg2namepd.set_index('id')['name'].to_dict()
#
# Reactions = []
# g = Graph()
# g.parse('rhea.rdf', format='xml')
# rh = Namespace("http://rdf.rhea-db.org/")
# rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
#
# # 查询反应和相关的化合物信息
# reaction_query = """
# SELECT DISTINCT ?x ?reaction ?general ?contain ?contains ?a ?n ?c ?f ?chebi ?cpd
# WHERE {
#     ?x rdfs:subClassOf rh:Reaction ;
#        rdfs:subClassOf ?general ;
#        rh:accession ?accession ;
#        rh:side ?reaction .
#
#     ?reaction ?contains ?contain .
#     ?contains rdfs:subPropertyOf rh:contains .
#     ?contain rh:compound ?cpd .
#     ?cpd rh:accession ?a .
#     }
# }
# """
# results = g.query(reaction_query)
#
# # 处理查询结果
# reaction_dict = {}
# for row in tqdm(results, desc="Processing Reactions"):
#     reaction = str(row.reaction)
#     if reaction not in reaction_dict:
#         reaction_dict[reaction] = {
#             'x': str(row.x),
#             'general': str(row.general),
#             'contains': []
#         }
#     if row.contain:
#         reaction_dict[reaction]['contains'].append((str(row.contain), str(row.contains)))
#
#     # 存储化合物信息
#     if row.a:
#         reaction_dict[reaction].setdefault('compounds', []).append({
#             'cpd': str(row.contain),
#             'accession': str(row.a),
#             'name': str(row.n),
#             'charge': int(row.c) if row.c else None,
#             'formula': str(row.f),
#             'chebi': str(row.chebi) if row.chebi else None
#         })
#
# # 处理包含 GENERIC 的化合物信息
# for row in tqdm(results, desc="Processing Reactions"):
#     if row.a and str(row.a).startswith('GENERIC'):
#         met = str(row.a)
#         cpdquery = f"""
#         SELECT DISTINCT ?cpd ?n ?c ?f ?chebi
#         WHERE {{
#             ?cpd rh:accession '{met}' ;
#                  rh:name ?n ;
#                  rh:reactivePart ?name .
#             ?name rh:charge ?c ;
#                   rh:formula ?f ;
#                   rh:chebi ?chebi.
#         }}
#         LIMIT 100
#         """
#         cpd_list = g.query(cpdquery)
#         for compound in cpd_list:
#             # 处理 GENERIC 化合物的相关信息
#             reaction_dict[str(row.reaction)].setdefault('compounds', []).append({
#                 'cpd': str(row.contain),
#                 'accession': str(compound.cpd),
#                 'name': str(compound.n),
#                 'charge': int(compound.c) if compound.c else None,
#                 'formula': str(compound.f),
#                 'chebi': str(compound.chebi) if compound.chebi else None
#             })
#
# # 添加反应到 cobra 模型
# for i, (x, data) in enumerate(reaction_dict.items()):
#     if i % 2 == 0:
#         reaction_add = cobra.Reaction()
#         result_dict = {}
#
#     for contain, contains in data['contains']:
#         for compound in data.get('compounds', []):
#             if contain == compound['cpd']:
#                 metabolite = cobra.Metabolite(
#                     name=compound['name'],
#                     charge=compound['charge'],
#                     id=compound['accession'],
#                     formula=compound['formula'],
#                     compartment='c'
#                 )
#                 if compound['chebi']:
#                     metabolite.annotation['chebi'] = compound['chebi'].split('_')[1]
#
#                 # 根据 contains 设置反应的方向
#                 if contains.split('/')[-1][8:] == 'N' or contains.split('/')[-1][8:] == 'n':
#                     result_dict[metabolite] = -4 if x.endswith('L') else 4
#                 elif contains.split('/')[-1][8:] == 'Nplus1':
#                     result_dict[metabolite] = -5 if x.endswith('L') else 5
#                 elif contains.split('/')[-1][8:] == 'Nminus1':
#                     result_dict[metabolite] = -3 if x.endswith('L') else 3
#                 elif contains.split('/')[-1][8:] == '2n':
#                     result_dict[metabolite] = -8 if x.endswith('L') else 8
#                 else:
#                     result_dict[metabolite] = -int(contains.split('/')[-1][8:]) if x.endswith('L') else int(contains.split('/')[-1][8:])
#
#     if i % 2 == 1:
#         reaction_add.lower_bound = 0
#         reaction_add.upper_bound = 1000
#         reaction_add.id = f'rhea_{data["x"].split("/")[-1]}'
#         reaction_add.annotation['rhea'] = data['x'].split("/")[-1]
#         reaction_add = add_anno(reaction_add, data['x'])
#         reaction_add.add_metabolites(result_dict)
#         Reactions.append(reaction_add)
# model=cobra.Model('universial')
# model.add_reactions(Reactions)
# cobra.io.save_json_model(model, 'reaction_pool.json')