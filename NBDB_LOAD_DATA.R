#Loads formatted data into Neo4j

library("RNeo4j")

#connect to Neo4j server
graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#add constrainst for each node
addConstraint(graph, "Agent", "name")
addConstraint(graph, "Pathway", "name")
addConstraint(graph, "CellLine", "name")
addConstraint(graph, "Combination", "name")
addConstraint(graph, "Gene", "name")

#add treatment data
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///TREATMENT_DATA.csv' AS row

MERGE (agent:Agent {name: UPPER(row.`Agent (name)`)}) 

MERGE (gene1:Gene {name: UPPER(row.`Gene1 (name)`)})

FOREACH(_ IN CASE WHEN row.`Gene2 (name)` <> '' THEN [1] ELSE [] END |
  MERGE (gene2:Gene {name: UPPER(row.`Gene2 (name)`)})
  MERGE (agent)-[:TARGETS]->(gene2)
)

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line (name)`)}) 

MERGE (agent)-[ic50:IC50]->(cellline)
SET ic50.score = toFloat(row.`IC50 (score)`),
ic50.test = row.`IC50 (test)`

MERGE (agent)-[:TARGETS]->(gene1)
;"
cypher(graph,query)

#add combination data
query = "
USING PERIODIC COMMIT //commits 1000 rows of imported data at a time
LOAD CSV WITH HEADERS FROM 'file:///COMBINATION_DATA.csv' AS row

//adds combination treatment nodes with name from file
MERGE (combination:Combination {name: [UPPER(row.`Agent1`), UPPER(row.`Agent2`)]})

//adds cell line data with name from file
MERGE (cellline:CellLine {name: UPPER(row.`Cell Line (name)`)})

//adds treatment agent with name from file
MERGE (agent1:Agent {name: UPPER(row.`Agent1`)})

//adds treatment agent with name from file
MERGE (agent2:Agent {name: UPPER(row.`Agent2`)})

//creates relationship from combination to single agent
MERGE (combination)-[:CONSISTS_OF]->(agent1)

//create relationship from combination to single agent
MERGE (combination)-[:CONSISTS_OF]->(agent2)

//create synergy relationships from combination to cell line with score, test, and strength as properties 
MERGE (combination)-[synergy:SYNERGY]->(cellline)
SET synergy.score = row.`SYNERGY (score)`,
synergy.test = toFloat(row.`SYNERGY (test)`),
synergy.strength = row.`SYNERGY (strength)`
;"
cypher(graph,query)

#add mutation data
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///MUT.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:MUTATED_IN]->(cellline)
;"
cypher(graph,query)

#add mutation annotation
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///MUT_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:MUTATED_IN]->(cellline)
;"
cypher(graph,query)

#add amplification annotation
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///AMP_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:AMPLIFIED_IN]->(cellline)
;"
cypher(graph,query)

#add deletion annotation
query = "
// adding deleted annotation
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///DEL_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:DELETED_IN]->(cellline)
;"
cypher(graph,query)

#add pathway data
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///PATH.csv' AS row

MERGE (pathway:Pathway {name: UPPER(row.`Pathway`)})

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (pathway)-[:ENCODED_BY]->(gene)
;"
cypher(graph,query)

#add expression outlier data
query = "
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///EXP_OUTLIERS.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)})

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:OVEREXPRESSED_IN]->(cellline)
;"
cypher(graph,query)
