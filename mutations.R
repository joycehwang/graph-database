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
