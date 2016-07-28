# Visualizations of subsets of the neutoblastoma graph database

library(RNeo4j)
library(igraph)
library(visNetwork)

#connect to Neo4j
graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#pathways associated with targeted gene
query = "
MATCH (agent:Agent)-[:TARGETS]->(gene:Gene)<-[:ENCODED_BY]-(pathway:Pathway)
RETURN agent.name AS from, pathway.name AS to, COUNT(*) AS weight
"
edges = cypher(graph,query) #edges data frame with from, to, and weight as variables

#nodes are "from" and "to" values from edges data frame
nodes = data.frame(id = unique(c(edges$from, edges$to)))
nodes$label = nodes$id

ig = graph_from_data_frame(edges, directed = FALSE)

nodes$value = betweenness(ig) #size of nodes function of betweenness centrality

clusters = cluster_edge_betweenness(ig) #colors nodes according to clusters
nodes$group = clusters$membership

visNetwork(nodes, edges)

#################################################################################

#cell lines tested by each agent
query = "
MATCH (agent:Agent)-[:IC50]->(cellline:CellLine)
RETURN agent.name AS from, cellline.name AS to, COUNT(*) AS weight
"
edges = cypher(graph,query) #edges data frame with from, to, and weight as variables

#nodes are "from" and "to" values from edges data frame
nodes = data.frame(id = unique(c(edges$from, edges$to)))
nodes$label = nodes$id

ig = graph_from_data_frame(edges, directed = FALSE)

nodes$value = betweenness(ig) #size of nodes function of betweenness centrality

clusters = cluster_edge_betweenness(ig) #colors nodes according to clusters
nodes$group = clusters$membership

visNetwork(nodes, edges)
