library(RNeo4j)

graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#provides list of relationships
summary(graph)

#opens up browser window in RStudio
browse(graph)

query = "
MATCH (mek1:Gene {name:'MAP2K1'})--(x)
RETURN x.name AS Node, labels(x) AS Type
;"
cypher(graph,query)

query = "
MATCH (mek2:Gene {name:'MAP2K2'})--(y)
RETURN y.name AS Node, labels(y) AS Type
;"
cypher(graph,query)

query = "
MATCH (erk:Pathway {name:'ST_ERK1_ERK2_MAPK_PATHWAY'})--(a)
RETURN a.name AS Node, labels(a) AS Type
;"
cypher(graph,query)

query = "
MATCH (mapk:Pathway {name:'KEGG_MAPK_SIGNALING_PATHWAY'})--(b)
RETURN b.name AS Node, labels(b) AS Type
;"
cypher(graph,query)

query = "
MATCH (rb1p52:Pathway {name:'MARTINEZ_RB1_AND_TP53_TARGETS_DN'})--(c)
RETURN c.name AS Node, labels(c) AS Type
;"
cypher(graph,query)

query = "
MATCH (p53:Pathway {name:'MARTINEZ_TP53_TARGETS_DN'})--(d)
RETURN d.name AS Node, labels(d) AS Type
;"
cypher(graph,query)

query = "
MATCH (rb1u:Pathway {name:'MARTINEZ_RB1_TARGETS_UP'})--(e)
RETURN e.name AS Node, labels(e) AS Type
;"
cypher(graph,query)

query = "
MATCH (rb1d:Pathway {name:'MARTINEZ_RB1_TARGETS_DN'})--(f)
RETURN f.name AS Node, labels(f) AS Type
;"
cypher(graph,query)

#######################################################################

#number of affected genes in each cell line
query = "
MATCH (g:Gene)-->(c:CellLine)
RETURN c.name AS CellLine, COUNT(distinct g) AS Gene
ORDER BY Gene DESC
"
data = cypher(graph,query)

#number of cell lines affected by each gene
query = "
MATCH (g:Gene)-->(c:CellLine)
RETURN g.name As Gene, COUNT(distinct c) AS CellLine
ORDER BY CellLine DESC
"
data = cypher(graph,query)
