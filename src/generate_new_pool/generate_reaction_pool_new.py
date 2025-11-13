from rdflib import Graph, Namespace, URIRef
import rdflib

def read_ttl_file(file_path: str):
    """
    读取并解析 .ttl 文件，打印 RDF 三元组
    :param file_path: .ttl 文件路径
    """
    file_path='./data_available/MNXmnet.ttl'
    g = Graph()
    try:
        # 解析 Turtle 文件
        g.parse(file_path, format="turtle")
        print(f"成功读取 {file_path}，共 {len(g)} 条三元组。\n")

    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 未找到。")
    except Exception as e:
        print(f"解析文件时出错: {e}")
    rm = Namespace(" https://rdf.metanetx.org/")
    rdfs = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    reaction_query = f"""
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX reac: <https://rdf.metanetx.org/reac/>
    PREFIX schema: <https://rdf.metanetx.org/schema/>
    SELECT ?rxns ?b
    WHERE {{
       ?rxns schema:mnxr 'https://rdf.metanetx.org/reac/MNXR174239' .
    }}
        """
    reactions = g.query(reaction_query)
    rxx=[]
    for row in reactions:
        rxx.append(row.b.split("/")[-1])
        print(row.b)

    rhea=[str(row.x).split("/")[-1] for row in reactions]


file_path='./data_available/MNXmnet.ttl'
g = Graph()
g.parse(file_path, format="turtle")