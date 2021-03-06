#CUSTOM SCRIPT FOR CREATING MATRIX OF ACCESSIBILITY PER CELL PER AGGREGATE PEAK
#version to execute from an IDE, not optimized for command line
import sys
import pysam
import numpy as np
import pandas as pd

numpeaks = 85044  

#open original sam file and barcode mapping file
samfile = open("SRR1947694_MQ_dedup_noM.sorted.sam","r")
barcodefile = open("GSM1647124_CtlSet1.indexconversion.txt","r")
newsamfile = open("SRR1947694_MQ_dedup_noM_barcodes_500.sorted.sam","w")

#create dictionary of barcodes: key is the read name, value is the "updated" barcode
#initialize dictionary of count of reads per cell: key is "updated" barcode, value is initialized to 0
barcodes ={};cellcount = {}
for line in barcodefile.readlines():
    line = line.split()
    barcodes["SRR1947694." + line[0]] = line[2]
    cellcount[line[2]] = 0
barcodefile.close()

missing = []
#loop through sam file and count how many reads correspond to each barcode, update dictionary
#start at line 87 to skip the header
for line in samfile.readlines()[87:]:
    line = line.split()
    try:
        cellcount[barcodes[line[0]]]+=1
    except:
        missing.append(line[0])

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
    try:
        if barcodes[line[0]] in keeplist:
            line[0] = barcodes[line[0]]
            line = "\t".join(line) + "\n"
            newsamfile.write(line)
    except:
        0

newsamfile.close()
samfile.close()

#have to stop here, add header, convert to bam, and make new index file

#open new sam file with pysam and peaks file
samfile = pysam.AlignmentFile("SRR1947694_MQ_dedup_noM_barcodes_500_header.sorted.bam", "rb")
peaks = open("SRR1947694_peaks.bed","r")

#intialize matrix to store data; x counts the number of cells
peakmatrix = np.zeros((len(keeplist), numpeaks+1),dtype=int)

#loop through peaks, store count in matrix
i = 0
chrs = []; starts = []; stops=[]
for peak in peaks.readlines():
    peak = peak.split()
    chrs.append(peak[0]); starts.append(peak[1]); stops.append(peak[2])
    #loop through all reads aligning to that peak
    for read in samfile.fetch(peak[0], int(peak[1]), int(peak[2])):
        cellnum = str(read).split()[0]
        #add a count to the matrix at the entry for cell # (row) and peak # (column), in the order they appear in the peak file
        peakmatrix[int(keeplist.index(cellnum)),i]+=1
    i+=1

peakmatrix = pd.DataFrame(peakmatrix,header=keeplist)
peakmatrix = peakmatrix.transpose()
peakmatrix['chr'] =  chrs; peakmatrix['start'] =  starts; peakmatrix['stop'] =  stops;
peakmatrix.to_csv("PeakMatrix_Ctl1.csv")
