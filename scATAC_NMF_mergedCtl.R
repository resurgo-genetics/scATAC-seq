library(NMF)
library(data.table)
library(proxy)
library(lsa)


setwd("~/Documents/LeslieLab/CusanovichData_Processed/CtlSetMerged_DHSmatrix_GM12878andHL60/")
#DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.csv',header=F,stringsAsFactors = F,sep =" ") 
#DHSmatrix = DHSmatrix[,4:ncol(DHSmatrix)]

#optional converstion to TF_IDF matrix
#DHSmatrix = lw_logtf(DHSmatrix) * gw_idf(DHSmatrix)

#on non-binary matrix
#NMF1 = nmf(DHSmatrix,4)
#saveRDS(NMF1,'CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4_tfidf.RDS')
NMF1 = readRDS('CtlSetMerged.dhsmatrix.fullyreduced.NMF2.RDS')
site_clusters = predict(NMF1,what='rows')
cell_clusters = predict(NMF1,what='columns')
lp<-coef(NMF1); mx<-basis(NMF1)
data<-mx%*%lp


#code to cluster data matrix with NMF-based clustering scheme
NMF2 = readRDS('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF7.RDS')
cell_clusters = predict(NMF2,what='columns')
data_clustered = cbind(site_clusters,data)
ordering1 = order(site_clusters)
mx = mx[ordering1,]
data_clustered = data_clustered[ordering1,]
data_clustered = rbind(cell_clusters,data_clustered[,2:ncol(data_clustered)])
ordering2 = order(cell_clusters)
lp = lp[,ordering2]
data_clustered = data_clustered[,ordering2]

write.table(data_clustered[2:nrow(data_clustered),],file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4_row7NMFclust.txt",sep="\t")
write.table(lp,"CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4.coeff.txt",row.names = F, quote=F,col.names = F)
write.table(mx,"CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4.basis.txt",row.names = F, quote=F,col.names = F)

#code to cluster data matrix with other methods
basisclust = hclust(dist(mx),method='ward.D')
coeffclust = hclust(dist(t(lp)),method='ward.D')
mx = mx[basisclust$order,]
data_clustered = data[basisclust$order,]
lp = lp[,coeffclust$order]
data_clustered = data_clustered[,coeffclust$order]

write.table(data_clustered[2:nrow(data_clustered),],file="CtlSetMerged.dhsmatrix.fullyreduced.NMF2hclust.txt",sep="\t")
write.table(lp,"CtlSetMerged.dhsmatrix.fullyreduced.NMF2hclust.coeff.txt",row.names = F, quote=F,col.names = F)
write.table(mx,"CtlSetMerged.dhsmatrix.fullyreduced.NMF2hclust.basis.txt",row.names = F, quote=F,col.names = F)


#code to calculate residuals based on objective function used in original algorithm
NMFresid = residuals(NMF1)

#code to sort cells by read depth
depth=read.delim('CtlSetMerged.dhsmatrix.fullyreduced.barcodes_readcounts.txt')
depthorder = order(depth$V2)
data_clustered = data[,depthorder]
lp = lp[,depthorder]

write.table(data_clustered,file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4_tfidf_readepth.txt",sep="\t")
write.table(lp,file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4_coeff_readepth.txt",sep="\t",col.names=F,row.names=F)

#code to sort cells by cell type assignment for mixture
cellorder = order(depth$V3)
data_clustered = data[,cellorder]
lp = lp[,cellorder]

write.table(data_clustered,file="CtlSetMerged.dhsmatrix.fullyreduced.NMF2_celltype.txt",sep="\t")
write.table(lp,file="CtlSetMerged.dhsmatrix.fullyreduced.NMF2_coeff_celltype.txt",sep="\t",col.names=F,row.names=F)

#sort the cells in order of the heirarchical clustering and print this out to a file
depth = depth[coeffclust$order,]
write.table(depth$V3,"CtlSetMerged.dhsmatrix.fullyreduced.NMF2_cellorderingfromhclust.txt",col.names=F,row.names=F,quote=F,na='0')
