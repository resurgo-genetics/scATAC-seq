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
ReadCounts = read.delim('SRR1947692_readcountpercell.txt',header=F,stringsAsFactors = F,col.names = c("Barcode","ReadCount"))
Assignments = merge(ReadCounts,Assignments,by='Barcode')
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]

#discrepancy between their read count and mine
plot(Assignments$ReadCount, Assignments$TotalReads, xlab = "My Read Count", ylab = "Their Read Count")

#calculate jaccard distances and scaling
distances = dist(PeakMatrix,method="jaccard",by_rows=T)
fit <- cmdscale(distances,eig=TRUE, k=2)

#check correlation (log transform)
x <- fit$points[,2] * -1
y = Assignments$ReadCount
plot(y,x,col = Assignments$CellTypeAssignment,xlab='log10(ReadCount))', ylab = 'Coordinate 1')
leg = unique(Assignments$CellTypeAssignment[!is.na(Assignments$CellTypeAssignment)])
legend('topleft', pch=c(1,1,1), col=1:3, legend=leg, bty='n', cex=.9, box.col='darkgreen')

#repeat for just GM12878 cells
PeakMatrix = PeakMatrix[(Assignments$CellTypeAssignment=="GM12878" & !is.na(Assignments$CellTypeAssignment)),]
Assignments = subset(Assignments, CellTypeAssignment=="GM12878")
distances = dist(PeakMatrix,method="jaccard",by_rows=T)
fit <- cmdscale(distances,eig=TRUE, k=10)

#check correlation (log transform)
x <- fit$points[,1] * -1
y = log10(Assignments$ReadCount)
plot(y,x,col = Assignments$CellTypeAssignment,xlab='log10(ReadCount))', ylab = 'Coordinate 1')
leg = unique(Assignments$CellTypeAssignment[!is.na(Assignments$CellTypeAssignment)])
legend('topleft', pch=c(1,1,1), col=1:3, legend=leg, bty='n', cex=.9, box.col='darkgreen')


#plot multidimensional scaling of subtypes
x = fit$points[,1]
y =  fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 subtypes, Jaccard, ATAC-seq Peaks")

#perform kmeans on scaled features
kmeans_GM = kmeans(fit$points[,1:10],centers=2,nstart=1)
clusters = kmeans_GM$cluster
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 subtypes, kmeans, ATAC-seq Peaks",col=clusters)
legend('topleft', pch=c(1,1,1), col=1:3, legend=c("Cluster 1","Cluster 2"), bty='n', cex=.9, box.col='darkgreen')

#perform kmeans on scaled features with C1 removed
kmeans_GM = kmeans(fit$points[,2:10],centers=2,nstart=1)
clusters = kmeans_GM$cluster
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 subtypes, kmeans, ATAC-seq Peaks",col=clusters)
legend('topleft', pch=c(1,1,1), col=1:3, legend=c("Cluster 1","Cluster 2"), bty='n', cex=.9, box.col='darkgreen')

#perform kmeans on all features
kmeans_GM = kmeans(PeakMatrix,centers=2,nstart=1)
clusters = kmeans_GM$cluster
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 subtypes, kmeans, ATAC-seq Peaks",col=clusters)
legend('topleft', pch=c(1,1,1), col=1:3, legend=c("Cluster 1","Cluster 2"), bty='n', cex=.9, box.col='darkgreen')
