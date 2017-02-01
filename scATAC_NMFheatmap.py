from matplotlib import pyplot as plt
import pylab as pl
import sys
import seaborn as sns
import numpy as np

inputfile = sys.argv[1]
#inputfile = "CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.NMF4clustered.txt"
outputfile = sys.argv[2]

NMF = open(inputfile,'r')
lines = NMF.readlines()
NMFmat = np.zeros([len(lines)-1,len(lines[0].split('\t'))])

i=0
for line in lines[1:]:
    NMFmat[i] = map(float,line.split()[1:])
    i+=1

#plt.clf()
#plt.ylabel('Sites')
#plt.xlabel('Cells')
fig = plt.imshow(NMFmat, cmap="seismic",vmin=0, vmax=.2,aspect=.5)
plt.axis('off')
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig(outputfile, bbox_inches='tight', pad_inches = 0,dpi=500)

#sns.clustermap(NMFmat,cmap='seismic',vmin=0, vmax=.2,aspect=.5)

