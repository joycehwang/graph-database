library(RNeo4j)

graph = startGraph("http://localhost:7474/db/data/", username = "neo4j", password = "NEO$J")

#number of BET sensitive cell lines (sentivity defined as IC50 < 500 nM)
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'})
WHERE r.score < 500
RETURN COUNT(c) AS count
"
SENS_CL = as.numeric(cypher(graph,query))

#number of BET resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'})
WHERE r.score > 500
RETURN COUNT(c) AS count
"
RES_CL = as.numeric(cypher(graph,query))

#################################################################################

#cell lines sensitive to BET inhibitors 
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'}), 
(g:Gene)-->(c)
WHERE r.score < 500
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
BET_SENS = cypher(graph,query)
BET_SENS$Count = BET_SENS$Count/SENS_CL

library(ggplot2)

ggplot(data = BET_SENS, aes(Count)) +
  geom_histogram(binwidth = 0.1) +
  geom_freqpoly(binwidth = 0.1) +
  labs(title = "Gene Count in BET Sensitive Cell Lines", x = "Gene Count", y = "Frequency")

#cell lines resistent to BET inhibitors 
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'}), 
(g:Gene)-->(c)
WHERE r.score > 500
RETURN g.name AS Gene, count(g) AS Count
ORDER BY Count DESC
;"
BET_RES = cypher(graph,query)
BET_RES$Count = BET_RES$Count/RES_CL

ggplot(data = BET_RES, aes(Count)) +
  geom_histogram(binwidth = 0.2) +
  geom_freqpoly(binwidth = 0.2) +
  labs(title = "Gene Count in BET Resistent Cell Lines", x = "Gene Count", y = "Frequency")

#################################################################################

#scatterplot of resistent counts vs sensitive counts for genes
intSets <- intersect(BET_RES[,1], BET_SENS[,1])
rownames(BET_RES) <- BET_RES[,1]
rownames(BET_SENS) <- BET_SENS[,1]
DF_RES_SENS_GENE <- cbind(BET_RES[intSets,], BET_SENS[intSets,2]);
colnames(DF_RES_SENS_GENE) <- c("Gene", "Resistant_Count", "Sensitive_Count");
grep(DF_RES_SENS_GENE[,1], "BRD");
DF_RES_SENS_GENE[,"BRD_SET"] <- grepl("BRD", DF_RES_SENS_GENE[,1]);
ggplot(DF_RES_SENS_GENE, aes(Resistant_Count, Sensitive_Count, color=BRD_SET)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw();
ggplot(DF_RES_SENS_GENE, aes(Resistant_Count, Sensitive_Count, color=BRD_SET, size=BRD_SET)) +
  geom_point() +
  theme_bw();

#sensitive:resistent counts ratio
DF_RES_SENS_GENE[,"ResSensRatio"] <- DF_RES_SENS_GENE["Resistant_Count"]/DF_RES_SENS_GENE[,"Sensitive_Count"];
DF_RES_SENS_GENE <- DF_RES_SENS_GENE[order(DF_RES_SENS_GENE[,"ResSensRatio"]),]

hist(DF_RES_SENS_GENE[,"ResSensRatio"], breaks = 10, main = "BET Resistant:Sensitive Gene Ratios", xlab = "Resistant:Sensitive Ratio")
DF_RES_SENS_GENE[DF_RES_SENS_GENE[,"BRD_SET"] == T,]

#log transform ratios to normalize data
LOG_TRANS_RATIO <- log(DF_RES_SENS_GENE$ResSensRatio)
BET_LOG_TRANS <- cbind(DF_RES_SENS_GENE$Gene, LOG_TRANS_RATIO)
hist(LOG_TRANS_RATIO, breaks = 10, main = "BET Resistant:Sensitive Gene Ratios", xlab = "ln(Resistant:Sensitive Ratio)")

#################################################################################

#pathways affected in sensitive cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'}),
(p:Pathway)-[*2]->(c)
WHERE r.score < 500
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
BET_SENS_PATH = cypher(graph,query)
BET_SENS_PATH$Count = BET_SENS_PATH$Count/SENS_CL

ggplot(data = BET_SENS_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in BET Sensitive Cell Lines", x = "Pathway Count", y = "Frequency")

#pathways affected in resistent cell lines
query = "
MATCH (c:CellLine)<-[r:IC50]-(a:Agent)-[:TARGETS]->(:Gene {name:'BRD4'}),
(p:Pathway)-[*2]->(c)
WHERE r.score > 500
RETURN p.name AS Pathway, count(p) AS Count
ORDER BY Count DESC
"
BET_RES_PATH = cypher(graph,query)
BET_RES_PATH$Count = BET_RES_PATH$Count/RES_CL

ggplot(data = BET_RES_PATH, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_freqpoly(binwidth = 1) +
  labs(title = "Pathway Count in BET Resistent Cell Lines", x = "Pathway Count", y = "Frequency")

#################################################################################

#scatterplot of resistent counts vs sensitive counts
intSets <- intersect(BET_RES_PATH[,1], BET_SENS_PATH[,1])
rownames(BET_RES_PATH) <- BET_RES_PATH[,1]
rownames(BET_SENS_PATH) <- BET_SENS_PATH[,1]
DF_RES_SENS_PATH <- cbind(BET_RES_PATH[intSets,], BET_SENS_PATH[intSets,2]);
colnames(DF_RES_SENS_PATH) <- c("Pathway", "Resistant_Count", "Sensitive_Count");
grep(DF_RES_SENS_PATH[,1], "BRD");
DF_RES_SENS_PATH[,"BRD_SET"] <- grepl("BRD", DF_RES_SENS_PATH[,1]);
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=BRD_SET)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title = "BET Sensitive and Resistant Pathway Counts", x = "Resistant Count", y = "Sensitive Count") +
  guides(color = guide_legend(title = "BET"));
ggplot(DF_RES_SENS_PATH, aes(Resistant_Count, Sensitive_Count, color=BRD_SET, size=BRD_SET)) +
  geom_point() +
  theme_bw() +
  labs(title = "BET Sensitive and Resistant Pathway Counts", x = "Resistant Count", y = "Sensitive Count") +
  guides(color = guide_legend(title = "BET"), size = guide_legend(title = "BET"));

#sensitive:resistent counts ratio
DF_RES_SENS_PATH[,"ResSensRatio"] <- DF_RES_SENS_PATH["Resistant_Count"]/DF_RES_SENS_PATH[,"Sensitive_Count"];
DF_RES_SENS_PATH <- DF_RES_SENS_PATH[order(DF_RES_SENS_PATH[,"ResSensRatio"]),]

hist(DF_RES_SENS_PATH[,"ResSensRatio"], breaks = 100, main = "BET Resistant:Sensitive Pathway Ratios", xlab = "Resistant:Sensitive Ratio")
DF_RES_SENS_PATH[DF_RES_SENS_PATH[,"BRD_SET"] == T,]

#################################################################################

#log transform ratios to normalize data
LOG_TRANS_RATIO <- log(DF_RES_SENS_PATH$ResSensRatio)
BET_LOG_TRANS <- cbind(DF_RES_SENS_PATH$Pathway, LOG_TRANS_RATIO)
hist(LOG_TRANS_RATIO, breaks = 150, main = "BET Resistant:Sensitive Pathway Ratios", xlab = "ln(Resistant:Sensitive Ratio)")