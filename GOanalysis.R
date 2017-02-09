#analysis and visualization of enriched GO terms over clusters
library(data.table)
library(wordcloud)
library(stringi)
GO1 = fread('cluster1GO.txt')
GO2 = fread('cluster2GO.txt')
GO3 = fread('cluster3GO.txt')
GO4 = fread('cluster4GO.txt')
GO5 = fread('cluster5GO.txt')

GO1_unique = subset(GO1,!(`GO biological process complete` %in% c(GO2$`GO biological process complete`,GO3$`GO biological process complete`,GO4$`GO biological process complete`,GO5$`GO biological process complete`)))
GO2_unique = subset(GO2,!(`GO biological process complete` %in% c(GO1$`GO biological process complete`,GO3$`GO biological process complete`,GO4$`GO biological process complete`,GO5$`GO biological process complete`)))
GO3_unique = subset(GO3,!(`GO biological process complete` %in% c(GO2$`GO biological process complete`,GO1$`GO biological process complete`,GO4$`GO biological process complete`,GO5$`GO biological process complete`)))
GO4_unique = subset(GO4,!(`GO biological process complete` %in% c(GO2$`GO biological process complete`,GO3$`GO biological process complete`,GO1$`GO biological process complete`,GO5$`GO biological process complete`)))
GO5_unique = subset(GO5,!(`GO biological process complete` %in% c(GO2$`GO biological process complete`,GO3$`GO biological process complete`,GO4$`GO biological process complete`,GO1$`GO biological process complete`)))

GO1_unique$`GO biological process complete` = stri_sub(GO1_unique$`GO biological process complete`,0,-14)
GO2_unique$`GO biological process complete` = stri_sub(GO2_unique$`GO biological process complete`,0,-14)
GO3_unique$`GO biological process complete` = stri_sub(GO3_unique$`GO biological process complete`,0,-14)
GO4_unique$`GO biological process complete` = stri_sub(GO4_unique$`GO biological process complete`,0,-14)
GO5_unique$`GO biological process complete` = stri_sub(GO5_unique$`GO biological process complete`,0,-14)


pdf("cluster1GO_uniqueterms.pdf",width=20,height=15)
wordcloud(GO1_unique$`GO biological process complete`,GO1_unique$`upload_1 (P-value)`,max.words =25,rot.per=0,color='red')
dev.off()

pdf("cluster2GO_uniqueterms.pdf",width=20,height=15)
wordcloud(GO2_unique$`GO biological process complete`,GO2_unique$`upload_1 (P-value)`,max.words =35,rot.per=0,color='blue')
dev.off()

pdf("cluster3GO_uniqueterms.pdf",width=20,height=15)
wordcloud(GO3_unique$`GO biological process complete`,GO3_unique$`upload_1 (P-value)`,max.words =35,rot.per=0,color='darkgreen')
dev.off()

pdf("cluster4GO_uniqueterms.pdf",width=20,height=15)
wordcloud(GO4_unique$`GO biological process complete`,GO4_unique$`upload_1 (P-value)`,max.words =35,rot.per=0,color='purple')
dev.off()

pdf("cluster5GO_uniqueterms.pdf",width=20,height=15)
wordcloud(GO5_unique$`GO biological process complete`,GO5_unique$`upload_1 (P-value)`,max.words =35,rot.per=0,color='orange')
dev.off()

#alternative form of visualization: heatmap of p-values
cluster_intersection=Reduce(intersect, list(GO1$`GO biological process complete`,GO2$`GO biological process complete`,GO3$`GO biological process complete`,GO4$`GO biological process complete`,GO5$`GO biological process complete`))
GO1=GO1[which(GO1$`GO biological process complete`%in%cluster_intersection),]
GO2=GO2[which(GO2$`GO biological process complete`%in%cluster_intersection),]
GO3=GO3[which(GO3$`GO biological process complete`%in%cluster_intersection),]
GO4=GO4[which(GO4$`GO biological process complete`%in%cluster_intersection),]
GO5=GO5[which(GO5$`GO biological process complete`%in%cluster_intersection),]


