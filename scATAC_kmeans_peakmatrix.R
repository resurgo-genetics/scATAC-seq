library(data.table)
library(proxy)
library(pheatmap)
library(RColorBrewer)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

PeakMatrix = read.delim('PeakMatrix_new.csv',header=T,stringsAsFactors = F,sep=',')
#have to edit the indexing on this once i have the data
barcodes = cbind(as.character(PeakMatrix[,1]),c(1:(nrow(PeakMatrix))))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
PeakMatrix = PeakMatrix[2:ncol(PeakMatrix)]
PeakMatrix[PeakMatrix>1]<-1

kmeans1 = kmeans(PeakMatrix,centers=2,nstart=1)
clusters = kmeans1$cluster

#plot kmeans clusters on a dimensionality reduction of the points in 2D space
distances = dist(PeakMatrix,method="jaccard",by_rows=T) #check if it actually is rows
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

Assignments = read.delim('GSM1647122_GM12878vsHEK.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 vs HEK293T, k-means, ATAC-seq Peaks",col=clusters)
legend('bottomleft', pch=c(1,1,1), col=1:3, c('Cluster 1', 'Cluster 2','Cluster 3'), bty='n', cex=.9, box.col='darkgreen')

#plot kmeans clusters together in a heatmap, cluster columns with heirarchical clustering based on jaccard distances
hc.cols <- hclust(dist(t(PeakMatrix),method='jaccard'))
PeakMatrixk = cbind(clusters,PeakMatrix)
ordering = order(PeakMatrixk[,1])
PeakMatrixk = PeakMatrix[ordering,]
PeakMatrix1 = subset(PeakMatrix,clusters==1)
PeakMatrix2 = subset(PeakMatrix,clusters==2)
png("test2.png",width=3000,height=1000)
heatmap(as.matrix(PeakMatrixk),Rowv=NA, Colv=as.dendrogram(hc.cols), col = heat.colors(256), scale="none",labRow = NA, labCol = NA,margins=c(1,.01))
dev.off()

#plot actual group assignments in a heatmap
PeakMatrixc = cbind(Assignments$CellTypeAssignment,PeakMatrix)
ordering = order(PeakMatrixc[,1])
PeakMatrixc = PeakMatrix[ordering,]
png("test2.png",width=3000,height=1000)
heatmap(as.matrix(PeakMatrixc),Rowv=NA, Colv=NA, col = heat.colors(256), scale="none",labRow = NA, labCol = NA,margins=c(1,.01))
dev.off()
