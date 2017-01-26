import numpy as np

inputfile = open("CtlSetMerged.dhsmatrix.txt",'r')
outputfile = open("CtlSetMerged.dhsmatrix.GM12878only_rowreduced.csv",'w')


columns = [0] * (len(inputfile.readline().split(','))-10)
inputfile.seek(0)
for line in inputfile.readlines()[1:]:
    line = np.array(line.split(','))
    try:
        line[np.where(line=="")] = '0'
        line = list(line[0:len(line)-1])
    except:
        line = list(line[0:len(line)-1])
    line1 = line[5:4542] + line[4544:9124] + line[9126:13841]
    line1 = map(int,line1)
    columns = np.add(columns,line1)
    if np.sum(line1) >= 150:
        outputfile.write("\t".join(line[1:4] + map(str,line1)) +'\n')
  
inputfile.close()
outputfile.close()

columns = np.array(columns)
columns = (columns >= 400)

inputfile = open("CtlSetMerged.dhsmatrix.GM12878only_rowreduced.csv",'r')
outputfile = open("CtlSetMerged.dhsmatrix.GM12878only_reduced.csv",'w')

for line in inputfile.readlines():
    line = map(int,line.split(',')))
    outputfile.write('\t'.join(line[columns]))
    
    
    

