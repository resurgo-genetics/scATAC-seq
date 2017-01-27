library(data.table)
library(proxy)
library(lsa)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

#original code
DHSmatrix = fread('CtlSetMerged.dhsmatrix.GMfullyreduced.csv',header=F,stringsAsFactors = F,sep ="\t")
barcodes = cbind(as.character(t(DHSmatrix[1,])),c(-4:(ncol(DHSmatrix)-5)))[6:ncol(DHSmatrix),]
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))
DHSmatrix = DHSmatrix[2:nrow(DHSmatrix),6:ncol(DHSmatrix)]
Assignments = read.delim('GSM1647124_CtlSet1.readcounts.txt',header=T)
Assignments = rbind(Assignments,read.delim('GSM1647125_CtlSet2.readcounts.txt',header=T))
Assignments = rbind(Assignments,read.delim('GSM1647126_CtlSet3.readcounts.txt',header=T))
Assignments = merge(barcodes,Assignments,by='Barcode',all.x=T)

#using LSA package; i don't think this mimics what is done in the paper exactly
tfmatrix = as.textmatrix(DHSmatrix)
weightedmatrix = lw_tf(tfmatrix) * gw_idf(tfmatrix)

#something to transform scores to -1.5 +1.5 scale?
weightedsvd = svd(weightedmatrix, nu = 6, nv = 6)


#my own code to perform LSA as closely as possible to what they do in the paper
  #create vector of total site usage per cell
  numsites = colSums(DHSmatrix)
  #divide each column by this value (or the log of it?)
  weightedDHSmatrix = DHSmatrix/colSums
  #weightedDHSmatrix = DHSmatrix/log10(colSums)
  
  


