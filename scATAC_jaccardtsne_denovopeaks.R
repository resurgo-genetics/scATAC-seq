library(data.table)
library(proxy)
library(Rtsne)
setwd("~/Documents/LeslieLab/CusanovichData_Processed/SRR1947692")

PeakMatrix = read.delim('PeakMatrix_SRR1947692.csv',header=T,stringsAsFactors = F,sep=',')
#have to edit the indexing on this once i have the data
barcodes = cbind(as.character(PeakMatrix[,1]),c(1:(nrow(PeakMatrix))))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
PeakMatrix = PeakMatrix[2:ncol(PeakMatrix)]
PeakMatrix[PeakMatrix>1]<-1

#may have to convert all to numerics, not sure yet  
distances = dist(PeakMatrix,method="jaccard",by_rows=T) #check if it actually is rows

tsne1 = Rtsne(distances, is_distance=T,verbose=T)

Assignments = read.delim('../GSM1647122_GM12878vsHEK.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]
legendnames = as.character(unique(Assignments$CellTypeAssignment))
legendnames = legendnames[!is.na(legendnames)]
plot(tsne1$Y[,1],tsne1$Y[,2], xlab="", ylab="",main="GM12878 vs. HEK293T, t-sne w/ Jaccard Distance",col=Assignments$CellTypeAssignment)
legend('bottomleft', pch=c(1,1,1), col=1:3, legendnames, bty='n', cex=.9, box.col='darkgreen')

