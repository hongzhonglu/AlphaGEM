#目前我们有ecnumber信息，我们需要对应到相应的反应，得到mapping关系1、gpr，2、reaction---comp，3、cpd。
from rdflib import Graph
g = Graph()

# 读取RDF文件
g.parse('rhea.rdf', format='xml')  # 假设RDF文件是XML格式

query="""
PREFIX rh: <http://rdf.rhea-db.org/> 
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT ?x ?y
WHERE{
?x rdfs:subClassOf rh:Reaction .
?x rh:accession "RHEA:10736" .
?x rdfs:label ?comment.
?x rh:bidirectionalReaction ?y .
}
"""

query2="""
PREFIX rh: <http://rdf.rhea-db.org/> 
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?reaction ?contain ?contains ?cpd3 ?a ?b ?n ?c
WHERE {
  {
    SELECT DISTINCT ?reaction
    WHERE {
      ?x rdfs:subClassOf rh:Reaction ;
         rh:ec <http://purl.uniprot.org/enzyme/3.6.1.9> ;
         rh:side ?reaction .
    }
    LIMIT 100
  }
  ?contains rdfs:subPropertyOf rh:contains .
  ?reaction ?contains ?contain .
  ?contain rh:compound ?cpd .
  ?cpd rh:chebi ?cpd3 ;
       rh:accession ?a ;
       rh:name ?n ;
       rh:charge ?b ;
       rh:formula ?c .
}
"""

results=g.query(query)

for row in results:
    print(row[0])