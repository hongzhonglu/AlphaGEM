from rdflib import Graph, Namespace, URIRef
import cobra
import pandas as pd


def check_reactions(reac,model):
    result={}
    for met in reac:
        metvalue=reac[met]
        metabolite=model.metabolites.get_by_id(met).elements
        result={k: result.get(k, 0) + metvalue*metabolite.get(k, 0) for k in set(result) | set(metabolite)}
    print(result)
    for i in result:
        if result[i]!=0:
            return 0
    return 1


def add_bireaction(g,enzyme_c,keggdict,model,gpr,model_met_exist,kegg2name,model_met_exist2):
    rh = Namespace("http://rdf.rhea-db.org/")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    reaction_query = f"""
    SELECT DISTINCT ?x ?reaction
    WHERE {{
      ?x rdfs:subClassOf rh:Reaction ;
         rh:ec <http://purl.uniprot.org/enzyme/{enzyme_c}> ;
         rh:bidirectionalReaction ?y ;
         rh:side ?reaction .
    }}
    LIMIT 100
    """
    reactions = g.query(reaction_query)

    # 将结果存储在一个字典中
    reaction_dict = {str(row.reaction): [] for row in reactions}
    if reaction_dict=={}:
        return model
    x_reaction_dict = {str(row.x): [] for row in reactions}
    for row in reactions:
        x_reaction_dict[str(row.x)].append(str(row.reaction))

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
    SELECT ?reaction ?contain ?contains ?cpd3 ?a ?b ?n ?c
    WHERE {
      VALUES ?contain { """ + " ".join(f"<{c[0]}>" for r in reaction_dict.values() for c in r) + """ }
      ?contain rh:compound ?cpd .
      ?cpd rh:chebi ?cpd3 ;
           rh:accession ?a ;
           rh:name ?n ;
           rh:charge ?b ;
           rh:formula ?c .
    }
    """
    compounds = g.query(compound_query)

    compound_info_dict = {}
    count = 0


    for x, reactions in x_reaction_dict.items():
        result_dict = {}
        reaction_add = cobra.Reaction()
        for reaction in reactions:
            contains_list = reaction_dict[reaction]
            for contain, contains in contains_list:
                for row in compounds:
                    if contain == str(row.contain):
                        metbolic = cobra.Metabolite(
                            name=f"{row.n}",
                            charge=int(f"{row.b}"),
                            id=f"{row.a}",
                            formula=f"{row.c}",
                            compartment='c'
                        )
                        metbolic.annotation['chebi'] = f"{row.a}"
                        try :
                            met=model_met_exist[f"{row.a}"]

                        except:
                            try:
                                id=model_met_exist2[f"{row.n}"]
                                if model.metabolites.get_by_id(id).compartment=='c' and model.metabolites.get_by_id(id).charge==int(f"{row.b}") and model.metabolites.get_by_id(id).formula==f"{row.c}":
                                    met=id
                                else:
                                    met = model_met_exist[f"{row.a}"]
                            except:
                                try:
                                    model.metabolites.get_by_id(f"{row.a}")
                                    met = f"{row.a}"
                                except:
                                    model.add_metabolites(metbolic)
                                    met = f"{row.a}"
                        if contains.split('/')[-1][8:] == 'N':
                            if reaction.endswith('L'):
                                result_dict[met] = -4
                            else:
                                result_dict[met] = 4
                        else:
                            if reaction.endswith('L'):
                                result_dict[met] = -int(contains.split('/')[-1][8:])
                            else:
                                result_dict[met] = int(contains.split('/')[-1][8:])
        reaction_add.lower_bound = 0
        reaction_add.upper_bound = 1000
        reaction_add.id = f'rhea_{x.split("/")[-1]}'
        reaction_add.annotation['rhea'] = x.split("/")[-1]
        try:
            reaction_add.annotation['kegg'] =keggdict[int(x.split("/")[-1])]
            try:
                reaction_add.name=kegg2name[reaction_add.annotation['kegg']].split(';')[0]
            except:
                1
        except:
            a=0
        if check_reactions(result_dict,model)==0:
            continue
        model.add_reactions([reaction_add])
        reaction_add.add_metabolites(result_dict)
        reaction_add.gene_reaction_rule=gpr
    return model


def model_reaction(name,refname):
    g = Graph()
    model = cobra.io.load_yaml_model(f'models/tarmodel{name}.yml')
    # 从文件或URL中解析RDF数据
    g.parse('rhea.rdf', format='xml')
    kegg2rhea = pd.read_csv('data_available/rhea2kegg_reaction.tsv', sep='\t')
    keggdict = kegg2rhea.set_index('RHEA_ID')['ID'].to_dict()
    model_met_exist = {}
    for met in model.metabolites:
        try:
            if type(met.annotation['chebi']) is list:
                for c in met.annotation['chebi']:
                    if met.compartment == 'c':
                        model_met_exist[c] = met.id
            else:
                if met.compartment == 'c':
                    model_met_exist[met.annotation['chebi']] = met.id
        except:
            1
    kegg2namepd = pd.read_csv('data_available/reaction', sep='\t', header=None, names=['id', 'name'])
    kegg2name = kegg2namepd.set_index('id')['name'].to_dict()
    if refname == 'yeast':
        model_met_exist['CHEBI:15378'] = 's_0794'
    model_met_exist2 = {}
    for met in model.metabolites:
        if met.compartment == 'c':
            model_met_exist2[met.name] = met.id
    gpr_e=pd.read_excel(f'juzhen/{name}gpregg_clean.xlsx')
    for index,row1 in gpr_e.iterrows():
        try:
            model = add_bireaction(g, row1['ec'], keggdict, model, row1['reaction'], model_met_exist, kegg2name, model_met_exist2)
        except:
            continue
    cobra.io.save_yaml_model(model, f'models/tarmodel_{name}.yml')
