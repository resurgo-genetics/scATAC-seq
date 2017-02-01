library(data.table)
library(proxy)
library(Rtsne)
setwd("~/Documents/LeslieLab/CusanovichData_Processed/CtlSetMerged_DHSmatrix_GM12878andHL60/")

#THIS FILE IS FOR CLUSTERING AND DIM. REDUCTION OF THE FULL HL-60 VS GM12878 SET USING ALL SITES
#get sites
DHSmatrix = fread('CtlSetMerged.dhsmatrix.fullyreduced.csv',header=F,stringsAsFactors = F,sep= " ")
sites = cbind(DHSmatrix[,1:3],c(1:nrow(DHSmatrix)))
setnames(sites,c('V1','V2','V3','order'))
annots = fread("CtlSetMerged.dhsmatrix.fullyreduced.annotated_sites.bed",header=F,stringsAsFactors = F)
annots = unique(annots,by=c('V1','V2','V3'))
annots = merge(annots,sites,by=c("V1",'V2','V3'))
annots = annots[order(as.numeric(annots$order)),]
rm(sites)

#get barcodes
DHSmatrix = fread('CtlSetMerged.dhsmatrix.rowreduced.csv',header=F,stringsAsFactors = F,sep ="\t")
DHSmatrix = DHSmatrix[1:nrow(DHSmatrix),4:ncol(DHSmatrix)]
keepcols = (colSums(DHSmatrix) >= 400)
barcodes = read.delim('CtlSetMerged.dhsmatrix.rowreduced.csv',nrows=1,header=F)[4:13836]
barcodes = cbind(t(barcodes),c(1:13833))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
barcodes = subset(barcodes,keepcols)
setnames(barcodes,c("Barcode","Order"))
Assignments = read.delim('../GSM1647124_CtlSet1.readcounts.txt',header=T)
Assignments = rbind(Assignments,read.delim('../GSM1647125_CtlSet2.readcounts.txt',header=T))
Assignments = rbind(Assignments,read.delim('../GSM1647126_CtlSet3.readcounts.txt',header=T))
Assignments = merge(barcodes,Assignments,by='Barcode',all.x=T)
Assignments = Assignments[order(as.numeric(Assignments$Order)),]
rm(DHSmatrix)
rm(barcodes)
rm(keepcols)

#t-sne on cells; replace distance file for non-HL-60 specific
distances = readRDS("CtlSetMerged.dhsmatrix.fullyreduced.cell_distances.RDS")
tsne1 = Rtsne(distances,is_distance=T,verbose=T)
plot(tsne1$Y[,1],tsne1$Y[,2], xlab = "",ylab='',main="Clustering of GM12878 and HL-60 Mixture by Site Usage",col=Assignments$CellTypeAssignment)

#t-sne on sites; replace distance file for non-HL-60 specific
distances2 = readRDS("")
tsne2 = Rtsne(distances2,is_distance=T,verbose=T)
plot(tsne2$Y[,1],tsne2$Y[,2], xlab = "",ylab='',col=annots$V8,main="Clustering of Accessible Sites by GM12878 Subtype Coordination")

#MDS on cells;
fit <- cmdscale(distances,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="Clustering of GM12878 Subtypes by Site Usage")

#MDS on sites
fit <- cmdscale(distances2,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="Clustering of Accessible Sites by GM12878 Subtype Coordination")


