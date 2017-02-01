library(NMF)
library(data.table)
library(proxy)
library(gplots)

#THIS SCRIPT IMPLEMENTS NMF ON GM12878 SUBTYPES 
setwd("~/Documents/LeslieLab/CusanovichData_Processed/CtlSetMerged_DHSmatrix_GM12878only/NoHLSites/")
#DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.csv',header=F,stringsAsFactors = F,sep =" ") 
#DHSmatrix = DHSmatrix[,4:ncol(DHSmatrix)]
#on non-binary matrix
#NMF1 = nmf(DHSmatrix,4)
#saveRDS(NMF1,'CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4.RDS')
NMF1 = readRDS('CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4.RDS')
site_clusters = predict(NMF1,what='rows')
cell_clusters = predict(NMF1,what='columns')
lp<-coef(NMF1); mx<-basis(NMF1)
data<-mx%*%lp

data_clustered = cbind(site_clusters,data)
ordering1 = order(site_clusters)
data_clustered = data_clustered[ordering1,]
data_clustered = rbind(cell_clusters,data_clustered[,2:ncol(data_clustered)])
ordering2 = order(cell_clusters)
data_clustered = data_clustered[,ordering2]

write.table(data_clustered[2:nrow(data_clustered),],file="CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4clustered.txt",sep="\t")
write.table(site_clusters[ordering1],"rowclusters.txt",row.names = F, quote=F,col.names = F)
write.table(cell_clusters[ordering2],"columnclusters.txt",row.names = F, quote=F,col.names = F)
