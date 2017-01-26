library(data.table)
library(proxy)


PeakMatrix = read.delim('PeakMatrix_3.csv',header=T,stringsAsFactors = F,sep=',')
#have to edit the indexing on this once i have the data
barcodes = cbind(as.character(PeakMatrix[,1]),c(1:(nrow(PeakMatrix))))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
PeakMatrix = PeakMatrix[2:ncol(PeakMatrix)]
PeakMatrix[PeakMatrix>1]<-1

#may have to convert all to numerics, not sure yet  
distances = dist(PeakMatrix,method="jaccard",by_rows=T) #check if it actually is rows
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

Assignments = read.delim('GSM1647123_GM12878vsHL.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]
legendnames = as.character(unique(Assignments$CellTypeAssignment))
legendnames = legendnames[!is.na(legendnames)]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 vs. HL-60, Jaccard Distance, ATAC-seq Peaks",col=Assignments$CellTypeAssignment)
legend('topright', pch=c(1,1), col=1:3, legendnames, bty='n', cex=.9, box.col='darkgreen')

