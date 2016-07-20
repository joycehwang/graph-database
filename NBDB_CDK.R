library(RNeo4j)

graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#number of CDK sensitive cell lines (sentivity defined as IC50 < 1000 nanomoles)
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'})
WHERE r.score < 1000
RETURN COUNT(c) AS count
"
SENS_CL = as.numeric(cypher(graph,query))

#number of CDK resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'})
WHERE r.score > 1000
RETURN COUNT(c) AS count
"
RES_CL = as.numeric(cypher(graph,query))

#cell lines sensitive to CDK inhibitors 
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'}), 
(g:Gene)-->(c)
WHERE r.score < 1000
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
CDK_SENS = cypher(graph,query)
CDK_SENS$Count = CDK_SENS$Count/SENS_CL

library(ggplot2)

ggplot(data = CDK_SENS, aes(Count)) +
  geom_histogram(binwidth = 0.1) +
  geom_freqpoly(binwidth = 0.1) +
  labs(title = "Gene Count in CDK Sensitive Cell Lines", x = "Gene Count", y = "Frequency")

#cell lines resistent to CDK inhibitors 
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'}), 
(g:Gene)-->(c)
WHERE r.score > 1000
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
CDK_RES = cypher(graph,query)
CDK_RES$Count = CDK_RES$Count/RES_CL

ggplot(data = CDK_RES, aes(Count)) +
  geom_histogram(binwidth = 0.2) +
  geom_freqpoly(binwidth = 0.2) +
  labs(title = "Gene Count in CDK Resistent Cell Lines", x = "Gene Count", y = "Frequency")

#pathways affected in sensitive cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'}),
(p:Pathway)-[*2]->(c)
WHERE r.score < 1000
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
CDK_SENS_PATH = cypher(graph,query)
CDK_SENS_PATH$Count = CDK_SENS_PATH$Count/SENS_CL

ggplot(data = CDK_SENS_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in CDK Sensitive Cell Lines", x = "Pathway Count", y = "Frequency")

#pathways affected in resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'CDK4'}),
(p:Pathway)-[*2]->(c)
WHERE r.score > 1000
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
CDK_RES_PATH = cypher(graph,query)
CDK_RENS_PATH$Count = CDK_RES_PATH$Count/RES_CL

ggplot(data = CDK_RES_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in CDK Resistent Cell Lines", x = "Pathway Count", y = "Frequency")

#scatterplot of resistent counts vs sensitive counts
intSets <- intersect(CDK_RES_PATH[,1], CDK_SENS_PATH[,1])
rownames(CDK_RES_PATH) <- CDK_RES_PATH[,1]
rownames(CDK_SENS_PATH) <- CDK_SENS_PATH[,1]
DF_RES_SENS_PATH <- cbind(CDK_RES_PATH[intSets,], CDK_SENS_PATH[intSets,2]);
colnames(DF_RES_SENS_PATH) <- c("Pathway", "Resistant_Count", "Sensitive_Count");
grep(DF_RES_SENS_PATH[,1], "CDK");
DF_RES_SENS_PATH[,"CDK_SET"] <- grepl("CDK", DF_RES_SENS_PATH[,1]);
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=CDK_SET)) +
    geom_point() +
    geom_smooth(method="lm") +
    theme_bw();
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=CDK_SET, size=CDK_SET)) +
    geom_point() +
    theme_bw();

#sensitive:resistent counts ratio
DF_RES_SENS_PATH[,"SensResRatio"] <- DF_RES_SENS_PATH["Resistant_Count"]/DF_RES_SENS_PATH[,"Sensitive_Count"];
DF_RES_SENS_PATH <- DF_RES_SENS_PATH[order(DF_RES_SENS_PATH[,"SensResRatio"]),]

hist(DF_RES_SENS_PATH[,"SensResRatio"], breaks=100)
DF_RES_SENS_PATH[DF_RES_SENS_PATH[,"CDK_SET"]==T,]
