library(lsa)
library(data.table)
library(proxy)

setwd("~/Documents/LeslieLab/CusanovichData_Processed/CtlSetMerged_DHSmatrix_GM12878only/NoHLSites/")
DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.csv',header=F,stringsAsFactors = F,sep =" ") 
DHSmatrix = DHSmatrix[,4:ncol(DHSmatrix)]

#optional converstion to TF_IDF matrix
#DHSmatrix = lw_logtf(DHSmatrix) * gw_idf(DHSmatrix)
#SVD1 = svd(DHSmatrix)


#saveRDS(SVD1, 'CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD.RDS')
SVD1 = readRDS('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD.RDS')
depth=6
u = SVD1$u; v = SVD1$v; d = diag(SVD1$d)
us <- as.matrix(u[, 2:depth])
vs <- as.matrix(v[, 2:depth])
ds <- as.matrix(d[2:depth,2:depth])
data <- us %*% ds %*% t(vs)

#code to cluster by ward heirarchical clustering
siteclust = hclust(dist(us),method='ward.D')
cellclust = hclust(dist(vs),method='ward.D')
data_clustered = data[siteclust$order,]
data_clustered = data_clustered[,cellclust$order]

write.table(data_clustered,file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD6_hclust.txt",sep="\t")

#code to sort sites by read depth
depth=read.delim('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.barcodes_readcounts.txt')
depthorder = order(depth$V2)
data_clustered = data[,depthorder]
vs = vs[depthorder,]

write.table(data_clustered,file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD6_hclust_readdepth.txt",sep="\t")
write.table(vs,file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD6_sortedv.txt",sep="\t",col.names=F,row.names=F)

#extract site clusters and write to bed file
png('test.png')
plot(siteclust)
dev.off()
siteclusters = cutree(siteclust,h=10)
sites = cbind(DHSmatrix[,1:3],rep("+",length(siteclusters)),siteclusters)
write.table(sites,"CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD2_6_siteclusters.bed",quote=F,row.names=F,sep='\t',col.names = F)

#calculate average within-cluster pairwise correlations
clustered_data = cbind(siteclusters,data)
cluster1_data = subset(clustered_data,siteclusters==1)
cluster2_data = subset(clustered_data,siteclusters==2)
cluster3_data = subset(clustered_data,siteclusters==3)
cluster4_data = subset(clustered_data,siteclusters==4)
cluster5_data = subset(clustered_data,siteclusters==5)

cluster1_cor = cor(t(cluster1_data[,2:ncol(cluster1_data)]))
cluster1_cor = cluster1_cor[lower.tri(cluster1_cor)]
cluster1_meancor = mean(cluster1_cor)

cluster2_cor = cor(t(cluster2_data[,2:ncol(cluster2_data)]))
cluster2_cor = cluster2_cor[lower.tri(cluster2_cor)]
cluster2_meancor = mean(cluster2x_cor)

cluster3_cor = cor(t(cluster3_data[,2:ncol(cluster3_data)]))
cluster3_cor = cluster3_cor[lower.tri(cluster3_cor)]
cluster3_meancor = mean(cluster3_cor)

cluster4_cor = cor(t(cluster4_data[,2:ncol(cluster4_data)]))
cluster4_cor = cluster4_cor[lower.tri(cluster4_cor)]
cluster4_meancor = mean(cluster4_cor)

cluster5_cor = cor(t(cluster5_data[,2:ncol(cluster5_data)]))
cluster5_cor = cluster5_cor[lower.tri(cluster5_cor)]
cluster5_meancor = mean(cluster5_cor)

all_cor = cor(t(data))
all_cor = all_cor[lower.tri(all_cor)]
all_cor = mean(all_cor)
