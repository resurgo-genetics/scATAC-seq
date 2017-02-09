from matplotlib import pyplot as plt
import pylab as pl
import sys
import seaborn as sns
import numpy as np

inputfile = sys.argv[1]
inputfile = "CtlSetMerged.dhsmatrix.GMfullyreduced_noHL.SVD6_hclust_readdepth.txt"
outputfile = sys.argv[2]

NMF = open(inputfile,'r')
lines = NMF.readlines()
NMFmat = np.zeros([len(lines)-1,len(lines[0].split())])

i=0
for line in lines[1:]:
    NMFmat[i] = map(float,line.split()[1:])
    i+=1

#NMFmat = (NMFmat-np.min(NMFmat))/(np.max(NMFmat)-np.min(NMFmat))

fig = plt.imshow(NMFmat, cmap="coolwarm",vmin=np.percentile(NMFmat,25),vmax = np.percentile(NMFmat,75),aspect=.5)
cbar = plt.colorbar(fig)
plt.axis('off')
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig(outputfile, bbox_inches='tight', pad_inches = 0,dpi=500)



