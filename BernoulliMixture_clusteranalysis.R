peaks = read.delim("SRR1947693_peaks.bed",header=F)
peaks = cbind(peaks,'index'=1:nrow(peaks))
peaks_annotated = read.delim("SRR1947693_peaks.annotated_sites.bed",header=F)
peaks_annotated = merge(peaks,peaks_annotated,by=c('V1','V2','V3'))

peaks_annotated$V8 = lapply(as.character(peaks_annotated$V8),function(x){
  x = strsplit(x,"_")
  x = x[[1]][length(x[[1]])]
  return(x)
})


#read in probabilities from each cluster and merge with peaks
cluster1 = read.delim('p1.txt',header=F)
cluster2 = read.delim('p2.txt',header=F)
clusters = cbind("Cluster1"=t(cluster1[1,]),"Cluster2"=t(cluster2[1,]))
clusters = cbind(clusters,"index"=1:nrow(clusters))
peaks_annotated_clusters = merge(peaks_annotated,clusters,by="index")

peaks_annotated_clusters_unique = unique(peaks_annotated_clusters,by="index")
#for each site, plot probability from each cluster against each other; color sites by type
colors = as.character(peaks_annotated_clusters$V8)
plot(as.numeric(peaks_annotated_clusters$`1`),as.numeric(peaks_annotated_clusters$`1.1`),col=as.factor(colors),pch=16,xlab="Bernoulli Parameters for Cell Type 1",ylab="Bernoulli Parameters for Cell Type 2")
legend('bottomright', pch=c(16), col=c(1:length(unique(colors))), unique(colors), bty='n', cex=.9, box.col='darkgreen')


#plot probabilities for each site as line graph, separate line for each cluster
#maybe sort by difference in probability
cluster1 = as.numeric(peaks_annotated_clusters$`1`)
cluster2 = as.numeric(peaks_annotated_clusters$`1.1`)
plot(1:length(cluster1),cluster1,cex=.000001,ylim=c(0,.4),col=2,xlab="")
lines(1:length(cluster1),cluster1,col=2)
par(new=T)
plot(1:length(cluster2),cluster2,cex=.000001,ylim=c(0,.4),col=3,xlab="Site",ylab="Bernoulli Parameter")
lines(1:length(cluster2),cluster2,col=3)
legend('topleft', pch=c(16), col=c(2,3), c("Cell Type 1", "Cell Type 2"), bty='n', cex=.9, box.col='darkgreen')
#maybe plot histogram of difference between probabilities

