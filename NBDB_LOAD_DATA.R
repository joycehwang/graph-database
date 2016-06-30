library("RNeo4j")

graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

addConstraint(graph, "Agent", "name")
addConstraint(graph, "Pathway", "name")
addConstraint(graph, "CellLine", "name")
addConstraint(graph, "Combination", "name")
addConstraint(graph, "Gene", "name")

query = "
// adding treatments
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
SET ic50.score = row.`IC50 (score)`,
ic50.test = row.`IC50 (test)`

MERGE (agent)-[:TARGETS]->(gene1)
;"
cypher(graph,query)

query = "
// adding combination treatments
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///COMBINATION_DATA.csv' AS row

MERGE (combination:Combination {name: [UPPER(row.`Agent1`), UPPER(row.`Agent2`)]})

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line (name)`)})

MERGE (agent1:Agent {name: UPPER(row.`Agent1`)})

MERGE (agent2:Agent {name: UPPER(row.`Agent2`)})

MERGE (combination)-[:CONSISTS_OF]->(agent1)

MERGE (combination)-[:CONSISTS_OF]->(agent2)

MERGE (combination)-[synergy:SYNERGY]->(cellline)
SET synergy.score = row.`SYNERGY (score)`,
synergy.strength = row.`SYNERGY (strength)`
;"
cypher(graph,query)

query = "
// adding mutations
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///MUT.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:MUTATED_IN]->(cellline)
;"
cypher(graph,query)

query = "
// adding mutations annotation
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///MUT_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:MUTATED_IN]->(cellline)
;"
cypher(graph,query)

query = "
// adding amplification annotation
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///AMP_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:AMPLIFIED_IN]->(cellline)
;"
cypher(graph,query)

query = "
// adding deleted annotation
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///DEL_ANN.csv' AS row

MERGE (cellline:CellLine {name: UPPER(row.`Cell Line`)}) 

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (gene)-[:DELETED_IN]->(cellline)
;"
cypher(graph,query)

query = "
// adding pathways
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM 'file:///PATH.csv' AS row

MERGE (pathway:Pathway {name: UPPER(row.`Pathway`)})

MERGE (gene:Gene {name: UPPER(row.`Gene`)})

MERGE (pathway)-[:ENCODED_BY]->(gene)
;"
cypher(graph,query)
