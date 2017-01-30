import numpy as np

inputfile = open("CtlSetMerged.dhsmatrix_dataonly.txt",'r')
outputfile = open("CtlSetMerged.dhsmatrix.rowreduced_gmonly.csv",'w')

GMsites = []
def subsetcells(line):
    out = []
    for i in line[5:]:
        if len(i) >=1 and len(i) < 5:
            out.append(int(i))
        elif len(i) == 0:
            out.append(0)
    return(out)
    
line = inputfile.readline().split(',')
columns = [0] * (len(line)-5); line1 = line[1:4]; line[-1] = line[-1][0:-1]
for i in line:
    if len(i) == 36:
        line1.append(i)
outputfile.write("\t".join(line1))


inputfile.seek(0)
for line in inputfile.readlines()[1:]:
    Gmbool=True
    line = line.split(',')
    line[-1] = line[-1][0:-1]
    try:
        line.index("Gmspecific")
    except:
        Gmbool=False
    line1 = subsetcells(line)
    columns = np.add(columns,np.array(line1))
    if np.sum(line1) >= 150 and Gmbool:
        outputfile.write("\t".join(line[1:4] + map(str,line1)) +'\n')
inputfile.close()
outputfile.close()

