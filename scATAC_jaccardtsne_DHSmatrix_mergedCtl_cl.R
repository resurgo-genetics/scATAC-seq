#!/usr/bin/env Rscript

library(data.table)
library(proxy)
library(Rtsne)

DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced.csv',header=F,stringsAsFactors = F,sep= " ")
sites = cbind(DHSmatrix$V1,DHSmatrix$V2,DHSmatrix$V3,c(1:nrow(DHSmatrix)))
sites = as.data.table(sites)
setnames(sites,c("V1",'V2','V3','order'))
DHSmatrix = subset(DHSmatrix, select = -c(V1,V2,V3))
DHSmatrix[DHSmatrix>1]<-1


annots = fread("CtlSetMerged.dhsmatrix.GMfullyreduced.annotated_sites.bed",header=F,stringsAsFactors = F)
annots = unique(annots,by=c('V1','V2','V3'))
sites$V2 = as.numeric(sites$V2); sites$V3 = as.numeric(sites$V3)
annots = merge(annots,sites,by=c("V1",'V2','V3'))
annots = annots[order(as.numeric(annots$order)),]

print("beginning distance calculation")
#on sites
distances = dist(DHSmatrix,method="jaccard",by_rows=F)
print("finished distance calculation; saving")
saveRDS("CtlSetMerged.dhsmatrix.GMfullyreduced.site_distances.RDS")
print("saved distance matrix to RDS file")

tsne1 = Rtsne(distances,is_distance=T,verbose=T)
pdf("GM12878subtypes_tnsesites.pdf", width=6, height=6)
plot(tsne1$Y[,1],tsne1$Y[,2], xlab = "",ylab='',main="Clustering of GM12878 Subtypes by Site Usage",col=annots$V8)
dev.off()
