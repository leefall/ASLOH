
par(mfrow=c(2,1),plt=c(0.05, 1, 0.3, 0.8))
output_corr_matrix <- read.delim("../../Data/TNBC_NGSCheckmate.txt")
data = output_corr_matrix
d3 <- as.dist((1 - data[,-1]))
clust3 <- hclust(d3, method = "average")
plot(clust3, lwd = 2, lty = 1,cex=0.45, xlab="Samples", sub = "",  ylab="Distance (1-Pearson correlation)",hang = -1, axes = FALSE,main="TNBC of TCGA")
axis(side = 2, at = seq(0, 1, 0.2), labels = FALSE, lwd = 2)
mtext(seq(0, 1, 0.2), side = 2, at = seq(0, 1, 0.2), line = 1,   las = 2)



output_corr_matrix <- read.delim("../../Data/mTNBC_NGSCheckmate.txt")
data = output_corr_matrix
d3 <- as.dist((1 - data[,-1]))
clust3 <- hclust(d3, method = "average")
plot(clust3, lwd = 2, lty = 1,cex=0.45, xlab="Samples", sub = "",  ylab="Distance (1-Pearson correlation)",hang = -1, axes = FALSE,main="mTNBC of Gustave Roussy")

