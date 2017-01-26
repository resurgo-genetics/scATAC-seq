#CUSTOM SCRIPT FOR CALCULATING PER-CELL READ COUNT FROM OUTPUT OF scATAC_peakmatrix_c1.py
#input file is a sam file with read names changed to cell barcodes of cells with at least 500 reads per cell
import sys

samfilename = sys.argv[1]

samfile = open(samfilename,"r")
outputfile = open("SRR1947693_readcountpercell.txt",'w')
cellcount = {}
for line in samfile.readlines()[87:]:
    line = line.split()
    try:
        cellcount[line[0]]+=1
    except:
        cellcount[line[0]]=1
samfile.close()
for code,count in cellcount.iteritems():
    line = code + '\t' + str(count) + '\n'
    outputfile.write(line)

outputfile.close()