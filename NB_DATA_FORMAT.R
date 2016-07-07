#reformatting cell line mutations data
df <- read.table("MUT_DATA_MATRIX.txt", header = TRUE)
col <- ncol(df)
cellLine <- colnames(df)
mut <- matrix()
mutations <- data.frame()
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  mut[i] <- subset(tmp, df[i] == 0.5)[1]
  if(length(mut[[i]]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(mut[[i]])), mut[i]))
    colnames(x) = NULL
    mutations <- as.matrix(rbind.data.frame(mutations, x))
  }
}
colnames(mutations) = c("Cell Line", "Gene")
write.csv(mutations, file = "MUT.csv", row.names = FALSE)

#reformatting cell line annotations data
df <- read.table("CELL_LINE_ANN.txt", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
cellLine <- colnames(df)
col <- ncol(df)
mut <- matrix()
mutations <- data.frame()
del <- matrix()
deleted <- data.frame()
amp <- matrix()
amplified <- data.frame()
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  mut[i] <- subset(tmp, df[i] == 0.5)[1]
  del[i] <- subset(tmp, df[i] == 1)[1]
  amp[i] <- subset(tmp, df[i] == 2)[1]
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
write.csv(mutations, file = "MUT_ANN.csv", row.names = FALSE)
write.csv(deleted, file = "DEL_ANN.csv", row.names = FALSE)
write.csv(amplified, file = "AMP_ANN.csv", row.names = FALSE)

#reformatting cell pathway data
col <- max(count.fields("c2.all.v5.1.symbols.gmt"))
df <- read.table("c2.all.v5.1.symbols.gmt", col.names = seq_len(col), na.strings=c("", "NA"), fill = TRUE)
df <- as.data.frame(t(df))
pathwayName <- df[1,]
col <- ncol(df)
pathways <- data.frame()
for (i in 1:col) {
  x = subset(df[i], !is.na(df[i]))[-c(1,2),]
  tmp <- as.matrix(cbind.data.frame(pathwayName[i], x))
  colnames(tmp) = NULL
  pathways <- as.matrix(rbind.data.frame(pathways, tmp))
}
colnames(pathways) = c("Pathway", "Gene")
write.csv(pathways, file = "PATH.csv", row.names = FALSE)

#reformatting expression outliers data
df <- read.table("ExpressionOutliersMatrix.txt", header = TRUE)
col <- ncol(df)
cellLine <- colnames(df)
out <- matrix()
outliers <- data.frame()
for (i in 2:col) {
  tmp <- data.frame(df[1], df[i])
  out[i] <- subset(tmp, df[i] == 1)[1]
  if(length(out[[i]]) > 0) {
    x = as.matrix(cbind.data.frame(rep(cellLine[i], times = length(out[[i]])), out[i]))
    colnames(x) = NULL
    outliers <- as.matrix(rbind.data.frame(outliers, x))
  }
}
colnames(outliers) = c("Cell Line", "Gene")
write.csv(outliers, file = "EXP_OUTLIERS.csv", row.names = FALSE)
