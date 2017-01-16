#CUSTOM SCRIPT FOR CREATING MATRIX OF ACCESSIBILITY PER CELL PER AGGREGATE PEAK
#user must input number of peaks in the peak file and cells as the first command line arguments
import sys
import pysam
import numpy as np
import pandas as pd

samname = sys.argv[1]
barcodename = sys.argv[2]
outputname = sys.argv[3]
SRAname = sys.argv[4]
keeplistfile = open("keepreads.txt",'w')

#open original sam file and barcode mapping file
samfile = open(samname,"r")
barcodefile = open(barcodename,"r")
newsamfile = open(outputname,"w")

#create dictionary of barcodes: key is the read name, value is the "updated" barcode
#initialize dictionary of count of reads per cell: key is "updated" barcode, value is initialized to 0
barcodes ={};cellcount = {}
for line in barcodefile.readlines():
    line = line.split()
    barcodes[SRAname + '.' + line[0]] = line[2]
    cellcount[line[2]] = 0
barcodefile.close()

#loop through sam file and count how many reads correspond to each barcode, update dictionary
#start at line 87 to skip the header
for line in samfile.readlines()[87:]:
    line = line.split()
    cellcount[barcodes[line[0]]]+=1

#determine which barcodes should be kept; cells with count > 500 and without ambiguous barcode sequences
#index of barcode in keeplist will serve as the row index in the accessibility matrix
keeplist = []
for code,count in cellcount.iteritems():
    if count > 500 and "_" not in code:
        keeplist.append(code); keeplistfile.write(code + '\n')
        
keeplistfile.close()

#to make more efficient later, just remove items from dictionary and use default function in following loop

samfile.seek(0)
#loop through sam file and replace readname with barcode, write to new sam file
for line in samfile.readlines()[87:]:
    line = line.split()
    if barcodes[line[0]] in keeplist:
        line[0] = barcodes[line[0]]
        line = "\t".join(line) + "\n"
        newsamfile.write(line)

newsamfile.close()
samfile.close()
