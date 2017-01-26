library(data.table)
library(proxy)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

#original code
#DHSmatrix = fread('CtlSetMerged.dhsmatrix.txt',header=F,stringsAsFactors = F,sep =",")
#saveRDS(DHSmatrix, 'CtlSetMerged.dhsmatrix.RDS')
#DHSmatrix = readRDS('CtlSetMerged.dhsmatrix.RDS')
DHSmatrix = readRDS('CtlSetMerged.dhsmatrix.GM12878only.RDS')
barcodes = cbind(as.character(t(DHSmatrix[1,])),c(-4:(ncol(DHSmatrix)-5)))[6:ncol(DHSmatrix),]
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
DHSmatrix = DHSmatrix[2:nrow(DHSmatrix),6:ncol(DHSmatrix)]
Assignments = read.delim('GSM1647124_CtlSet1.readcounts.txt',header=T)
Assignments = rbind(Assignments,read.delim('GSM1647125_CtlSet2.readcounts.txt',header=T))
Assignments = rbind(Assignments,read.delim('GSM1647126_CtlSet3.readcounts.txt',header=T))
Assignments = merge(barcodes,Assignments,by='Barcode',all.x=T)


#remove cells with < 400 and sites with < 150
#if you want to subset just cell types
#DHSmatrix = subset(DHSmatrix,select=(Assignments$CellTypeAssignment=="GM12878" & !is.na(Assignments$CellTypeAssignment)))
#rm(DHSmatrix)
#saveRDS(DHSmatrix1, 'CtlSetMerged.dhsmatrix.GM12878only.RDS')
DHSmatrix1 = sapply(DHSmatrix, as.numeric)
DHSmatrix1 = DHSmatrix[rowSums(DHSmatrix)>=150,]
DHSmatrix2 = DHSmatrix1[,colSums(DHSmatrix1)>=400]
saveRDS(DHSmatrix2, 'CtlSetMerged,dhsmatrix.GM12878only.reduced.RDS')
distances = dist(DHSmatrix2,method="jaccard",by_rows=F)
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 Subtypes, Jaccard Distance, DHS Peaks")
legend('bottomleft', pch=c(1), col=1, 'GM12878', bty='n', cex=.9, box.col='darkgreen')
