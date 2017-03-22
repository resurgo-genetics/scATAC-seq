import numpy as np
from matplotlib import pyplot as plt
import sklearn.feature_extraction
import pandas as pd
import seaborn
from scipy.sparse import csr_matrix
import lda

peakmatrix = pd.read_csv("PeakMatrix_SRR1947695.csv",index_col=0)
peakmatrix = peakmatrix.drop("#pattern",axis=1)
#map motifs to their transcription factors
TFdict2 = {}
lines = open("TF_Information_all_motifs_plus.txt",'r').readlines()
for line in lines[0:-2]:
    line = line.split()
    TFdict2[line[3]] = line[6]
TFindex = [TFdict2[motifid] for motifid in list(peakmatrix)]
peakmatrix.columns = TFindex

peakmatrix = peakmatrix.groupby(lambda x:x, axis=1).sum()

tfidf= sklearn.feature_extraction.text.TfidfTransformer()
peakmatrix_normalized = pd.DataFrame(tfidf.fit_transform(peakmatrix).todense(),columns=TFindex)


#visualization of cell by transcription factor data
    #number of cells each TF is present in
peakmatrix_binary = sklearn.preprocessing.binarize(peakmatrix)
TFsums = np.sum(peakmatrix_binary,axis=0)
plt.hist(TFsums,bins=100)

#remove "stop words"
stopwords = np.array(TFindex)[np.where(TFsums>=3500)]
peakmatrix_reduced = peakmatrix.drop(stopwords,axis=1)

model = lda.LDA(n_topics = 2, n_iter=1500)
model.fit(np.array(peakmatrix_reduced))

topic_word = model.topic_word_
n = 20
for i, topic_dist in enumerate(topic_word):
    topic_words = np.array(list(peakmatrix))[np.argsort(topic_dist)][:-(n+1):-1]
    print('*Topic {}\n- {}'.format(i, '\n'.join(topic_words)))


#clustering of cells
    #get cell types
celltypes = {}
barcodefile = open("/Users/cassandraburdziak/Documents/scATACseq//GSE67446_RAW/GSM1647125_CtlSet2.readcounts.txt","r").readlines()
for line in barcodefile[1:]:
    line = line.split()
    celltypes[line[0]] = line[6]
Cellindex = np.array([celltypes[barcode] for barcode in peakmatrix.index])


doc_topic = model.doc_topic_
cols = np.array(Cellindex)
cols[np.where(Cellindex=="GM12878")] = .01
cols[np.where(Cellindex=="HL60")] = .99
cols[np.where(Cellindex=="Mixed")] = .5
cols[np.where(Cellindex=="NA")] = .5
#cols = pd.DataFrame(cols,index=peakmatrix.index)
#doc_topic = pd.DataFrame(doc_topic,index=peakmatrix.index)
seaborn.clustermap(doc_topic,row_colors=cols,xticklabels=False,)

#visualization of results
    #what are the most prominent transcription factors in each topic: stacked bar chart

    #what are the most prominent topics in each cluster: stacked bar chart
    
    #topic by TF heatmap
    
    #pathway analysis of terms in each topic
    
    #do the clusters separate cell type: heatmap with annotation for both cell type and cluster membership
    
    