library(data.table)
library(proxy)
library(Rtsne)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced.csv',header=F,stringsAsFactors = F,sep= " ")
sites = cbind(DHSmatrix[,1:3],c(1:nrow(DHSmatrix)))
setnames(sites,c("V1",'V2','V3','order'))
DHSmatrix = DHSmatrix[,4:ncol(DHSmatrix)]
DHSmatrix[DHSmatrix>1]<-1

annots = fread("CtlSetMerged.dhsmatrix.GMfullyreduced.annotated_sites.bed",header=F,stringsAsFactors = F)
annots = unique(annots,by=c('V1','V2','V3'))
annots = merge(annots,sites,by=c("V1",'V2','V3'))
annots = annots[order(as.numeric(annots$order)),]

#on cells
distances = dist(DHSmatrix,method="jaccard",by_rows=F)
tsne1 = Rtsne(distances,is_distance=T,verbose=T)
plot(tsne1$Y[,1],tsne1$Y[,2], xlab = "",ylab='',main="Clustering of GM12878 Subtypes by Site Usage")

#on sites
distances = dist(DHSmatrix,method="jaccard",by_rows=T)
tsne2 = Rtsne(distances,is_distance=T,verbose=T)
plot(tsne2$Y[,1],tsne2$Y[,2], xlab = "",ylab='',col=annots$V8,main="Clustering of Accessible Sites by GM12878 Subtype Coordination")

