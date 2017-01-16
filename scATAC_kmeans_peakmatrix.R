library(data.table)
library(proxy)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

PeakMatrix = read.delim('PeakMatrix_new.csv',header=T,stringsAsFactors = F,sep=',')
#have to edit the indexing on this once i have the data
barcodes = cbind(as.character(PeakMatrix[,1]),c(1:(nrow(PeakMatrix))))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
PeakMatrix = PeakMatrix[2:ncol(PeakMatrix)]
PeakMatrix[PeakMatrix>1]<-1

kmeans1 = kmeans(PeakMatrix,centers=3,nstart=1)
clusters = kmeans1$cluster

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
