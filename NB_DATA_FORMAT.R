# Format raw data into file with just the direct relationships to read into Neo4j

#################################################################################

#reformatting cell line mutations data
#raw data in the form of matrix with cell lines and genes as the axes
df <- read.table("MUT_DATA_MATRIX.txt", header = TRUE) 
col <- ncol(df)
#cell line names
cellLine <- colnames(df)
mut <- matrix()
mutations <- data.frame()
#genes with mutations for each cell line
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  #values of 0.5 indicates mutation
  mut[i] <- subset(tmp, df[i] == 0.5)[1]
  #if there are mutations present add set of gene and cell line relationships to data frame
  if(length(mut[[i]]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(mut[[i]])), mut[i]))
    colnames(x) = NULL
    mutations <- as.matrix(rbind.data.frame(mutations, x))
  }
}
colnames(mutations) = c("Cell Line", "Gene")
#export as csv
write.csv(mutations, file = "MUT.csv", row.names = FALSE)

#################################################################################

#reformatting cell line annotations data
#raw data in the form of matrix with cell lines and genes as the axes
df <- read.table("CELL_LINE_ANN.txt", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
#cell line names
cellLine <- colnames(df)
col <- ncol(df)
mut <- matrix()
mutations <- data.frame()
del <- matrix()
deleted <- data.frame()
amp <- matrix()
amplified <- data.frame()
#genes with mutations, amplifications, deletions for each cell line
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  #values of 0.5 indicate mutation
  mut[i] <- subset(tmp, df[i] == 0.5)[1]
  #values of 1 indicate deletion
  del[i] <- subset(tmp, df[i] == 1)[1]
  #values of 2 indicate amplification
  amp[i] <- subset(tmp, df[i] == 2)[1]
  #if any modifications are present add gene and cell line relationship to data frame
  if (length(mut[[i]]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(mut[[i]])), mut[i]))
    colnames(x) = NULL
    mutations <- as.matrix(rbind.data.frame(mutations, x))
  } 
  if (length(del[i]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(del[[i]])), del[i]))
    colnames(x) = NULL
    deleted <- as.matrix(rbind.data.frame(deleted, x))
  } 
  if (length(amp[i]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(amp[[i]])), amp[i]))
    colnames(x) = NULL
    amplified <- as.matrix(rbind.data.frame(amplified, x))
  } 
}
colnames(mutations) = c("Cell Line", "Gene")
colnames(deleted) = c("Cell Line", "Gene")
colnames(amplified) = c("Cell Line", "Gene")
#export as csv
write.csv(mutations, file = "MUT_ANN.csv", row.names = FALSE)
write.csv(deleted, file = "DEL_ANN.csv", row.names = FALSE)
write.csv(amplified, file = "AMP_ANN.csv", row.names = FALSE)

#################################################################################

#reformatting cell pathway data
#raw data as pathway name, url, list of associated genes
col <- max(count.fields("c2.all.v5.1.symbols.gmt"))
df <- read.table("c2.all.v5.1.symbols.gmt", col.names = seq_len(col), na.strings=c("", "NA"), fill = TRUE)
df <- as.data.frame(t(df))
#pathway names
pathwayName <- df[1,]
col <- ncol(df)
pathways <- data.frame()
for (i in 1:col) {
  #subset genes in each pathway
  x = subset(df[i], !is.na(df[i]))[-c(1,2),]
  #combine gene and pathway relationship
  tmp <- as.matrix(cbind.data.frame(pathwayName[i], x))
  colnames(tmp) = NULL
  pathways <- as.matrix(rbind.data.frame(pathways, tmp))
}
colnames(pathways) = c("Pathway", "Gene")
#export as csv
write.csv(pathways, file = "PATH.csv", row.names = FALSE)

#################################################################################

#reformatting expression outliers data
#raw data in the form of matrix with cell lines and genes as the axes
df <- read.table("ExpressionOutliersMatrix.txt", header = TRUE)
col <- ncol(df)
#cell line names
cellLine <- colnames(df)
out <- matrix()
outliers <- data.frame()
#genes with expression outliers for each cell lines
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  out[i] <- subset(tmp, df[i] == 1)[1]
  #if outliers are present add gene and cell line relationship to data frame
  if(length(out[[i]]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(out[[i]])), out[i]))
    colnames(x) = NULL
    outliers <- as.matrix(rbind.data.frame(outliers, x))
  }
}
colnames(outliers) = c("Cell Line", "Gene")
#export as csv
write.csv(outliers, file = "EXP_OUTLIERS.csv", row.names = FALSE)
