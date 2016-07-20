library(RNeo4j)

graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#number of MEK sensitive cell lines (IC50 < 1.16 micromole)
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'})
WHERE r.score < 1160
RETURN COUNT(c) AS count
"
SENS_CL = as.numeric(cypher(graph,query))

#number of MEK resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'})
WHERE r.score > 1160
RETURN COUNT(c) AS count
"
RES_CL = as.numeric(cypher(graph,query))

#cell lines sensitive to MEK inhibitors 
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'}), 
(g:Gene)-->(c)
WHERE r.score < 1160
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
MEK_SENS = cypher(graph,query)
MEK_SENS$Count = MEK_SENS$Count/SENS_CL

library(ggplot2)

ggplot(data = MEK_SENS, aes(Count)) +
  geom_histogram(binwidth = 0.1) +
  geom_freqpoly(binwidth = 0.1) +
  labs(title = "Gene Count in Sensitive Cell Lines", x = "Gene Count", y = "Frequency")
mean <- mean(MEK_SENS$Count)
sd <- sd(MEK_SENS$Count)

#cell lines resistent to MEK inhibitors
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'}),
(g:Gene)-->(c)
WHERE r.score > 1160
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
MEK_RES = cypher(graph,query)
MEK_RES$Count = MEK_RES$Count/RES_CL

ggplot(data = MEK_RES, aes(Count)) +
  geom_histogram(binwidth = 0.1) +
  geom_freqpoly(binwidth = 0.1) +
  labs(title = "Gene Count in Resistent Cell Lines", x = "Gene Count", y = "Frequency")
mean <- mean(count)
sd <- sd(count)

#pathways affected in sensitive cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'}),
(p:Pathway)-[*2]->(c)
WHERE r.score < 1160
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
MEK_SENS_PATH = cypher(graph,query)
MEK_SENS_PATH$Count = MEK_SENS_PATH$Count/SENS_CL

ggplot(data = MEK_SENS_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in Sensitive Cell Lines", x = "Pathway Count", y = "Frequency")

#pathways affected in resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'MAP2K1'}),
(p:Pathway)-[*2]->(c)
WHERE r.score > 1160
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
MEK_RES_PATH = cypher(graph,query)
MEK_RES_PATH$Count = MEK_RES_PATH$Count/RES_CL

ggplot(data = MEK_RES_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in Resistent Cell Lines", x = "Pathway Count", y = "Frequency")

#scatterplot of resistent counts vs sensitive counts
intSets <- intersect(MEK_RES_PATH[,1], MEK_SENS_PATH[,1])
rownames(MEK_RES_PATH) <- MEK_RES_PATH[,1]
rownames(MEK_SENS_PATH) <- MEK_SENS_PATH[,1]
DF_RES_SENS_PATH <- cbind(MEK_RES_PATH[intSets,], MEK_SENS_PATH[intSets,2]);
colnames(DF_RES_SENS_PATH) <- c("Pathway", "Resistant_Count", "Sensitive_Count");
grep(DF_RES_SENS_PATH[,1], "MAPK");
DF_RES_SENS_PATH[,"MAPK_SET"] <- grepl("MAPK", DF_RES_SENS_PATH[,1]);
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=MAPK_SET))+geom_point() +
  geom_smooth(method="lm") +
  theme_bw();
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=MAPK_SET, size=MAPK_SET)) +
  geom_point() +
  theme_bw();

#sensitive:resistent counts ratio
DF_RES_SENS_PATH[,"SensResRatio"] <- DF_RES_SENS_PATH["Resistant_Count"]/DF_RES_SENS_PATH[,"Sensitive_Count"];
DF_RES_SENS_PATH <- DF_RES_SENS_PATH[order(DF_RES_SENS_PATH[,"SensResRatio"]),]

hist(DF_RES_SENS_PATH[,"SensResRatio"], breaks=1000)
DF_RES_SENS_PATH[DF_RES_SENS_PATH[,"MAPK_SET"]==T,]

#look at BIOCARTA_RACCYCD_PATHWAY
query = "
MATCH (p:Pathway {name:'BIOCARTA_RACCYCD_PATHWAY'})-->(g:Gene)-[r]->(c:CellLine)
RETURN g.name AS Gene, type(r) AS Relationship, c.name AS CellLine
"
RACCYCD = cypher(graph,query)
