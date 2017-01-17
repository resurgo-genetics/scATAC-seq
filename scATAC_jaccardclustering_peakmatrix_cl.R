#!/usr/bin/env Rscript

library(data.table)
library(proxy)

args = commandArgs(trailingOnly=TRUE)
peakmatrixname = args[1]
readcountsname = args[2]
figuretitle = args[3]

PeakMatrix = read.delim(peakmatrixname,header=T,stringsAsFactors = F,sep=',')
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

Assignments = read.delim(readcountsname,header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
Assignments$Order = as.numeric(Assignments$Order)
Assignments = Assignments[order(Assignments$Order),]
legendnames = as.character(unique(Assignments$CellTypeAssignment))
legendnames = legendnames[!is.na(legendnames)]
pdf(figuretitle + ".pdf",width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main=figuretitle,col=Assignments$CellTypeAssignment)
legend('bottomleft', pch=c(2,2), col=1:3, legendnames, bty='n', cex=.9, box.col='darkgreen')
dev.off()
