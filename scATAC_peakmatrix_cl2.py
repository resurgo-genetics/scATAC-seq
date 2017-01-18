#CUSTOM SCRIPT FOR CREATING MATRIX OF ACCESSIBILITY PER CELL PER AGGREGATE PEAK
#user must input number of peaks in the peak file and cells as the first command line arguments
import sys
import pysam
import numpy as np
import pandas as pd
import getopt

args = sys.argv[1:]
bamname,peaksname,numpeaks = getopt.getopt(args,"b:p:n:")[0]
bamname,peaksname,numpeaks = bamname[1],peaksname[1],numpeaks[1]

samfile = pysam.AlignmentFile(bamname, "rb")
peaks = open(peaksname,"r")
keeplistfile = open("truecells.txt",'r')

numcells = len(keeplistfile.readlines())
keeplistfile.seek(0); keeplist = []
for i in keeplistfile.readlines():
    keeplist.append(i.split()[0])

#intialize matrix to store data; x counts the number of cells
peakmatrix = np.zeros((int(numcells), int(numpeaks)+1),dtype=int)

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
