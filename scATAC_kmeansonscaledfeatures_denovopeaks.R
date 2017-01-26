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
Assignments = read.delim('GSM1647122_GM12878vsHEK.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]

#perform dimensionality reduction on the data to 100 dimensions
distances = dist(PeakMatrix,method="jaccard",by_rows=T) #check if it actually is rows
fit <- cmdscale(distances,eig=TRUE, k=10)
reducedPeakMatrix <- fit$points[,1:10]
reducedPeakMatrix_2 <- fit$points[,2:10]

kmeans1 = kmeans(reducedPeakMatrix,centers=2,nstart=1)
clusters1 = kmeans1$cluster

kmeans2 = kmeans(reducedPeakMatrix_2,centers=2,nstart=1)
clusters2 = kmeans2$cluster

#plot kmeans clusters together in a heatmap
hc.cols <- hclust(dist(t(reducedPeakMatrix),method='jaccard'))
k = cbind(clusters1,reducedPeakMatrix)
ordering = order(k[,1])
k = reducedPeakMatrix[ordering,]
heatmap(as.matrix(k),Rowv=NA, Colv=NA, col = heat.colors(256), scale="row",labRow = NA, labCol = NA,margins=c(1,.01))

#plot actual group assignments in a heatmap
c = cbind(Assignments$CellTypeAssignment,reducedPeakMatrix)
ordering = order(c[,1])
c = reducedPeakMatrix[ordering,]
heatmap(as.matrix(c),Rowv=NA, Colv=as.dendrogram(hc.cols), col = heat.colors(256), scale="row",labRow = NA, labCol = NA,margins=c(1,.01))

#plot kmeans clusters together in a heatmap
hc.cols <- hclust(dist(t(reducedPeakMatrix_2),method='jaccard'))
k = cbind(clusters2,reducedPeakMatrix)
ordering = order(k[,1])
k = reducedPeakMatrix_2[ordering,]
heatmap(as.matrix(k),Rowv=NA, Colv=NA, col = heat.colors(256), scale="row",labRow = NA, labCol = NA,margins=c(1,.01))

#plot actual group assignments in a heatmap
c = cbind(Assignments$CellTypeAssignment,reducedPeakMatrix)
ordering = order(c[,1])
c = reducedPeakMatrix[ordering,]
heatmap(as.matrix(c),Rowv=NA, Colv=NA, col = heat.colors(256), scale="row",labRow = NA, labCol = NA,margins=c(1,.01))

#plot on the first two coordinates
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 vs HEK293T, kmeans, Coordinate 1 Included",col=clusters1)
legend('bottomleft', pch=c(1,1,1), col=1:3, c('Cluster 1', 'Cluster 2'), bty='n', cex=.9, box.col='darkgreen')

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 vs HEK293T, kmeans, Coordinate 1 Not Included",col=clusters2)
legend('bottomleft', pch=c(1,1,1), col=1:3, c('Cluster 1', 'Cluster 2'), bty='n', cex=.9, box.col='darkgreen')
