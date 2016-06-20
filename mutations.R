df <- read.table("MUT_DATA_MATRIX.txt", header = TRUE)
col <- ncol(df)
row <- nrow(df)
cellLine <- colnames(df)
mutations <- matrix()
for (i in 2:col) {
  tmp <- data.frame(df[1],df[i])
  mutations[i] <- subset(tmp,df[i] == 0.5)[1]
  write.table(mutations[i], file = paste(cellLine[i], ".csv", sep = ""), row.names = FALSE, col.names = FALSE)
}
