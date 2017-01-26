library(logisticPCA)
library(ggplot2)
library(data.table)
library(tsne)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")
DHSmatrix = read.delim('GSM1647121_GM12878vsPatski.dhsmatrix.txt',header=F,stringsAsFactors = F)
barcodes = cbind(as.character(t(DHSmatrix[1,])),c(-3:(ncol(DHSmatrix)-4)))[5:537,]
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
Assignments = read.delim('GSM1647121_GM12878vsPatski.readcounts.txt',header=T)
Assignments = merge(barcodes,Assignments,by='Barcode')
DHSmatrix = DHSmatrix[2:nrow(DHSmatrix),5:ncol(DHSmatrix)]
DHSmatrix = sapply(DHSmatrix, as.numeric)

#traditional PCA
pca1 = prcomp(t(DHSmatrix),scale=T,center=T)
scores = data.frame(pca1$x[,1:3])
plot(scores$PC1,scores$PC2,col = Assignments$SpeciesAssignment,main="PCA of GM12878 vs. Patski", xlab="PC1",ylab="PC2")

#logistic PCA; takes too long
pca = logisticPCA(DHSmatrix,k=2,m=4,quiet=F)
plot(logisticPCA,type = "scores") + geom_point(aes(colour = Assignments$CellTypeAssignment))

#t-sne
tsne1 = tsne(DHSmatrix)
