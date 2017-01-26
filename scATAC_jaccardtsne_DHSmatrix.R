library(data.table)
library(proxy)
library(Rtsne)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

#original code
DHSmatrix = read.delim('GSM1647122_GM12878vsHEK.dhsmatrix.txt',header=F,stringsAsFactors = F)
barcodes = cbind(as.character(t(DHSmatrix[1,])),c(-3:(ncol(DHSmatrix)-4)))[5:704,]
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
DHSmatrix = DHSmatrix[2:nrow(DHSmatrix),5:ncol(DHSmatrix)]
DHSmatrix[DHSmatrix>1]<-1
Assignments = read.delim('GSM1647122_GM12878vsHEK.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
legendnames = as.character(unique(Assignments$CellTypeAssignment))
legendnames = legendnames[!is.na(legendnames)]

DHSmatrix1 = sapply(DHSmatrix, as.numeric)
distances = dist(DHSmatrix,method="jaccard",by_rows=F)

#on distances
tsne1 = Rtsne(distances,is_distance=T,verbose=T)
plot(tsne1$Y[,1],tsne1$Y[,2], xlab = "",ylab='',main="Clustering of Accessible Sites")
legend('bottomleft', pch=c(1,1,1),legend=legendnames, col=1:3, , bty='n', cex=.9, box.col='darkgreen')

