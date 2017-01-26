import pandas
strings = ['chr','start','end','annot']

one = open("GSM1647124_CtlSet1.dhsmatrix.txt","r")
barcodesone = {key: object for key in one.readline().split()}
one.close()
two = open("GSM1647125_CtlSet2.dhsmatrix.txt","r")
barcodestwo = {key: object for key in two.readline().split()}
two.close()
three = open("GSM1647126_CtlSet3.dhsmatrix.txt","r")
barcodesthree =  {key: object for key in three.readline().split()}
three.close()

for dict in [barcodesone,barcodestwo,barcodesthree]:
    for key in strings:
        dict[key] = 'str'


one = pandas.read_csv("GSM1647124_CtlSet1.dhsmatrix.txt",sep="\t",dtype=barcodesone)
two = pandas.read_csv("GSM1647125_CtlSet2.dhsmatrix.txt",sep="\t",dtype=barcodestwo)
print "loaded data"
onetwo = pandas.merge(one,two,on=["chr","start","end"],how='outer')
onetwo = onetwo.fillna(value=0)
del one
del two
print "merged database"

onetwo.to_csv("CtlSetMerged_half.dhsmatrix.txt",index_label=False,header=False,sep='\t')
print "wrote one two to file"
three = pandas.read_csv("GSM1647126_CtlSet3.dhsmatrix.txt",sep="\t",dtype=barcodesthree)

onetwothree = pandas.merge(onetwo,three,on=["chr","start","end"],how='outer')
onetwothree = onetwothree.fillna(value=0)
onetwothree.to_csv("CtlSetMerged.dhsmatrix.txt",index_label=False,header=False)