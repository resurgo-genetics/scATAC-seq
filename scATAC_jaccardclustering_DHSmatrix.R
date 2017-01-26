library(data.table)
library(proxy)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

#original code
DHSmatrix = read.delim('GSM1647123_GM12878vsHL.dhsmatrix.txt',header=F,stringsAsFactors = F)
barcodes = cbind(as.character(t(DHSmatrix[1,])),c(-3:(ncol(DHSmatrix)-4)))[5:704,]
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
DHSmatrix = DHSmatrix[2:nrow(DHSmatrix),5:ncol(DHSmatrix)]
DHSmatrix[DHSmatrix>1]<-1


DHSmatrix1 = sapply(DHSmatrix, as.numeric)
distances = dist(DHSmatrix1,method="jaccard",by_rows=F)
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

Assignments = read.delim('GSM1647123_GM12878vsHL.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')


plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 vs HL-60, Manhattan Distance",col=Assignments$CellTypeAssignment,ylim = c(-1000,300), xlim =c(-150,1000))
legend('bottomleft', pch=c(2,2), col=1:3, c('GM12878', 'HL-60','Mixed'), bty='n', cex=.9, box.col='darkgreen')

#if you want to subset just cell types
#have to fix this, it is elminating rows instead of columns
#is this also a problem with the other script i have for my own peaks?
DHSmatrix2 = DHSmatrix[,(Assignments$CellTypeAssignment=="GM12878" & !is.na(Assignments$CellTypeAssignment))]
DHSmatrix2 = sapply(DHSmatrix2, as.numeric)
distances = dist(DHSmatrix2,method="jaccard",by_rows=F)
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="GM12878 Subtypes, Jaccard Distance, DHS Peaks")
legend('bottomleft', pch=c(1), col=1, 'GM12878', bty='n', cex=.9, box.col='darkgreen')

#plotting correlation between read depth and coordinate 1
y = -log10(Assignments$TotalReads)
plot(x2,y,xlab="log10(Read Depth)", ylab="Coordinate 1")
