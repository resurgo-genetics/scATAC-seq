library(data.table)
setwd("~/Documents/LeslieLab/CusanovichData_Processed")

DHSmatrix = fread('CtlSetMerged.dhsmatrix.rowreduced.csv',header=F,stringsAsFactors = F,sep ="\t")
DHSmatrix = DHSmatrix[1:nrow(DHSmatrix),4:ncol(DHSmatrix)]
barcodes = read.delim('CtlSetMerged.dhsmatrix.rowreduced_fixed.csv',nrows=1,header=F)[4:13836]
barcodes = cbind(t(barcodes),c(1:13833))
barcodes = as.data.frame(barcodes,stringsAsFactors = F)
setnames(barcodes,c("Barcode","Order"))


#subsetting to remove anything non-GM12878
Assignments = read.delim('GSM1647124_CtlSet1.readcounts.txt',header=T)
Assignments = rbind(Assignments,read.delim('GSM1647125_CtlSet2.readcounts.txt',header=T))
Assignments = rbind(Assignments,read.delim('GSM1647126_CtlSet3.readcounts.txt',header=T))
Assignments = merge(barcodes,Assignments,by='Barcode',all.x=T)
Assignments = Assignments[order(as.numeric(Assignments$Order)),]
#DHSmatrix = subset(DHSmatrix,select=(Assignments$CellTypeAssignment=="GM12878" & !is.na(Assignments$CellTypeAssignment)))
keepcols = (colSums(DHSmatrix) >= 400)
DHSmatrix = subset(DHSmatrix,select=keepcols)

#write to a file
write.table(DHSmatrix,"CtlSetMerged.dhsmatrix.fullyreduced.csv")

#get bed file for coordinates
DHSmatrix1 = fread('CtlSetMerged.dhsmatrix.rowreduced.csv',header=F,stringsAsFactors = F,sep ="\t",skip=1)
DHSmatrix1 = cbind(DHSmatrix1[,1:3],rep("+", 35526))
write.table(DHSmatrix1,"CtlSetMerged.dhsmatrix.GMfullyreduced.sites.bed",quote=F,row.names = F,col.names = F,sep="\t")

#write to file with site coordinates included
DHSmatrix = cbind(DHSmatrix1[,1:3],DHSmatrix)
write.table(DHSmatrix,"CtlSetMerged.dhsmatrix.GMfullyreduced.csv",row.names=F,quote=F,col.names=F)
