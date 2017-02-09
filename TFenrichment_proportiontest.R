library(data.table)
library(stringi)
library(gplots)
setwd("~/Documents/LeslieLab/CusanovichData_Processed/CtlSetMerged_DHSmatrix_GM12878only/NoHLSites/")

#read in total sites bed file to get total number of sites in each cluster
metadata = fread("/Users/cass/Documents/LeslieLab/hg19_annotations/ChIP_GM12878metadata.tsv")
clusters.bed = fread("CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD2_6_siteclusters.bed")
cluster.counts = table(clusters.bed$V5)
rm(clusters.bed)
#read in bed file with transcription factor annotations
clusters.TF.bed = fread("CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD2_6_siteclusters_CHiP.bed")
#subset by overlap, choose > 50 bp
clusters.TF.bed = subset(clusters.TF.bed, V17 >=50)

#get list of the names of TF file names for only hg19 experiments
clusters.TF.bed$'File accession' = stri_sub(clusters.TF.bed$V6,0,-8)
clusters.TF.bed = merge(clusters.TF.bed,metadata,by='File accession')
TF.names = unique(subset(clusters.TF.bed,Assembly=='hg19')$'Experiment target')

#define function to perform Fisher's exact test with a table of counts of sites per cluster
#the full bed with TF names, and the name of TF of interest
#output will be a list of p-values for enrichment in each cluster
TFEnrichment = function(TF.name,cluster.counts,TF.df){
  TF.df = subset(TF.df,`Experiment target`==TF.name)
  TF.df = unique(TF.df,by=c('V1','V2','V3'))
  clusters = names(cluster.counts)
  cluster.TF.counts = table(TF.df$V5)
  pvals = rep(0,length(clusters))
  for(i in 1:length(clusters))
  {PropCluster = cluster.TF.counts[i]/cluster.counts[i]
  PropNonCluster = sum(cluster.TF.counts)/sum(cluster.counts)
  pvals[i] = t.test(c(PropCluster,PropNonCluster),alternative="greater",paired=F)$p.value
  }
  return(pvals)}


#calculate p-values for each
pvalues = matrix(0, ncol = 6, nrow = length(TF.names))
pvalues = data.frame(pvalues)
for(x in 1:length(TF.names)){
  pvalues[x,] = c(TF.names[x],TFEnrichment(TF.names[x],cluster.counts,clusters.TF.bed))
}

#bonferroni correction
pvalues_correct = sapply(pvalues[,2:6],as.numeric)*(5*length(TF.names))

#benjamini hochberg correction instead
FDR = .05
pvals_sorted = sort(as.vector(t(pvalues[,2:6])))
limits = FDR*as.vector(1:length(pvals_sorted))/length(pvals_sorted)
significant_pvalues = pvalues[,2:6] < pvals_sorted[48]
qvals = matrix(p.adjust(as.matrix(pvalues[,2:6]),'fdr'),ncol=5,byrow=T)

#find significant enrichments
significant_pvalues = (pvalues_correct <= .05)

#get them
cluster1 = pvalues$X1[which(significant_pvalues[,1]==T)]
cluster2 = pvalues$X1[which(significant_pvalues[,2]==T)]
cluster3 = pvalues$X1[which(significant_pvalues[,3]==T)]
cluster4 = pvalues$X1[which(significant_pvalues[,4]==T)]
cluster5 = pvalues$X1[which(significant_pvalues[,5]==T)]

cluster1=stri_sub(cluster1,0,-7)
cluster2 = stri_sub(cluster2,0,-7)
cluster3 = stri_sub(cluster3,0,-7)
cluster4 =stri_sub(cluster4,0,-7)
cluster5 =stri_sub(cluster5,0,-7)

write.table(cluster1,"cluster1_enrichedTFs_BH.txt",quote=F,col.names=F,row.names=F)
write.table(cluster2,"cluster2_enrichedTFs_BH.txt",quote=F,col.names=F,row.names=F)
write.table(cluster3,"cluster3_enrichedTFs_BH.txt",quote=F,col.names=F,row.names=F)
write.table(cluster4,"cluster4_enrichedTFs_BH.txt",quote=F,col.names=F,row.names=F)
write.table(cluster5,"cluster5_enrichedTFs_BH.txt",quote=F,col.names=F,row.names=F)

#sanity check for tail of hypothesis test
MYCtot = length(which(clusters.TF.bed$`Experiment target`=="MYC-human"))
MYCtot/nrow(clusters.TF.bed)

oneMYC = length(which(subset(clusters.TF.bed,`Experiment target`=='MYC-human')$V5==1))
onetot = length(which(clusters.TF.bed$V5==1))
oneMYC/onetot #larger

twoMYC = length(which(subset(clusters.TF.bed,`Experiment target`=='MYC-human')$V5==2))
twotot = length(which(clusters.TF.bed$V5==2))
twoMYC/twotot #smaller


#heatmap of TF's by q-values
pdf('TFenrichmentheatmap_BH.pdf' ,width = 4,height=6)
heatmap.2(qvals,margins = c(5,5),labRow=stri_sub(TF.names,0,-7),density.info='none',key.xlab='q-value',cexRow = .2,cexCol=.5,labCol=NA,trace='none',dendrogram='row',key.title='',keysize = 1.5)
dev.off()


