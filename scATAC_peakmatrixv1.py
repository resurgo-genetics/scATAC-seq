#CUSTOM SCRIPT FOR CREATING MATRIX OF ACCESSIBILITY PER CELL PER AGGREGATE PEAK
#version to execute from an IDE, not optimized for command line
import sys
import pysam
import numpy as np
import pandas as pd

numpeaks = 38995 

#open original sam file and barcode mapping file
samfile = open("SRR1947692_MQ_dedup_noM.sorted.sam","r")
barcodefile = open("GSM1647122_GM12878vsHEK.indexconversion.txt","r")
newsamfile = open("SRR1947692_MQ_barcodes_dedup_500.sorted.sam","w")

#create dictionary of barcodes: key is the read name, value is the "updated" barcode
#initialize dictionary of count of reads per cell: key is "updated" barcode, value is initialized to 0
barcodes ={};cellcount = {}
for line in barcodefile.readlines():
    line = line.split()
    barcodes["SRR1947692." + line[0]] = line[2]
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
        keeplist.append(code)

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

#have to stop here, add header, convert to bam, and make new index file

#open new sam file with pysam and peaks file
samfile = pysam.AlignmentFile("SRR1947692_MQ_barcodes_dedup_noM_500_header.sorted.bam", "rb")
peaks = open("SRR1947692_peaks_new.bed","r")

#intialize matrix to store data; x counts the number of cells
peakmatrix = np.zeros((len(keeplist), numpeaks+1),dtype=int)

#loop through peaks, store count in matrix
i = 0
for peak in peaks.readlines():
    peak = peak.split()
    #loop through all reads aligning to that peak
    for read in samfile.fetch(peak[0], int(peak[1]), int(peak[2])):
        cellnum = str(read).split()[0]
        #add a count to the matrix at the entry for cell # (row) and peak # (column), in the order they appear in the peak file
        peakmatrix[int(keeplist.index(cellnum)),i]+=1
    i+=1


df = pd.DataFrame(peakmatrix,index=keeplist)
df.to_csv("PeakMatrix_new.csv")
