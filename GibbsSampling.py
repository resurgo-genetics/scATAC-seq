#COLLAPSED GIBBS SAMPLER USING CRP OF DIRICHLET PROCESS MIXTURE MODEL
import math
import numpy as np
import scipy.misc
import pickle
from matplotlib import pyplot as plt

def save_obj(obj, name):
    with open('obj\ '+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def load_obj(name ):
    with open('obj\ ' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

#read in data
    #store as numpy array where each row contains the data for one cell
data = open("/Users/cass/Documents/LeslieLab/CusanovichData_Processed/SRR1947693/PeakMatrix_SRR1947693.csv","r")
lines=data.readlines()

data=np.zeros(shape=(len(lines),len(lines[1].split(',')[1:])))
i=0
for line in lines[1:]:
    data[i] = map(int,line.split(",")[1:]); i+=1

#train_csv = open("/Users/cass/Documents/AppliedMachineLearning/HW1/train.csv","r")
#lines = train_csv.readlines()[1:]
#data = np.zeros([len(lines),len(lines[0].split(','))-1])

#truth = []; i=0
#for line in lines:
#    data[i] = map(float,line.split(',')[1:])
#    truth.append(line[0])
#    i+=1
#truth = map(int,truth)

#just take a subset
#data = data[0:1000]
#truth = truth[0:1000]

data[np.where(data>=1)]=1

#Gibbs Sampler Functions
#function to update parameter values
def update_bernoullip(zs,x):
    pnew = {}
    #for a current value of z, update parameters for each cell in the data matrix
    for j in set(zs):
        jcells=np.where(zs==j)[0]
        pnew[j] = np.sum(x[jcells],axis=0)/float(len(jcells))
    return pnew
        
#function to sample from posterior distribution of z given conditional probabilities
#idea from http://www.hongliangjie.com/2016/03/12/tricks-in-sampling-discrete-distributions-1-sampling-logspace/
#fix tomorrow
def sample_z(zvals, z_conditionals):
    weights = [z_conditionals[0]]
    for i in range(1,len(z_conditionals)):
        weights.append(weights[i-1]+z_conditionals[i])
    U = [np.random.uniform(0,weights[-1])]*len(weights)
    U = np.subtract(weights,U)
    return zvals[np.where(U==min(i for i in U if i > 0))[0][0]]

def sample_z_log(zvals, z_conditionals):
    weights = [z_conditionals[0]]
    for i in range(1,len(z_conditionals)):
        weights.append(np.logaddexp(weights[i-1],z_conditionals[i]))
    U = np.log(np.random.uniform(0,1)) + weights[-1]
    for k in range(len(zvals)):
        if U < weights[k]:
            break
    return zvals[k]

#function to calculate conditional probabilities of z
def conditprob_zj(zs_minusn,x_minusn,x_n,hyperbeta,hypergamma,hyperalpha):
    z_conditionals = []; zvals=[]; N=x_minusn.shape[0];
    #probabilities for existing tables
    for j in set(zs_minusn):
        jcells=np.where(zs_minusn==j)[0]; x_minusn_j = x_minusn[jcells]
        N_minusn_j = x_minusn_j.shape[0]; N_j = len(jcells)
        probx1 = 0; probx2 = 0; y4=math.log(hyperbeta/(hyperbeta+hypergamma)); y5=math.log(hypergamma/(hyperbeta+hypergamma))
        for d in range(x_n.shape[0]):
            x = x_n[d]
            y1=np.sum(x_minusn_j[:,d])
            y2=math.pow((hyperbeta+y1),x)
            y3=math.pow(hypergamma+N_j-y1,1-x)
            probx1=probx1+math.log(y2*y3/(hyperbeta+hypergamma+N_j))
            probx2= probx2+(y4*x+y5*(1-x))
        probz = math.log(N_minusn_j / (float(N) - 1 + hyperalpha))
        z_conditionals.append(probz+probx1); zvals.append(j)
    probz = math.log(hyperalpha/(N-1+hyperalpha))
    z_conditionals.append(probz+probx2); zvals.append(j + 1)
    return(z_conditionals,zvals)
    


#random initialization of z_j's and p_j's 
z = np.random.choice([1,2],size=data.shape[0])

#initialize hyperparameters
hyperbeta=.5;hypergamma=.5;hyperalpha = 100.0

#for each iteration
for i in range(10):
    #for each data point
    print i
    for n in range(data.shape[0]):
        #unassign data point
        x_n = data[n]
        x_minusn = np.concatenate((data[0:n],data[n+1:]),axis=0)
        zs_minusn = np.concatenate((z[0:n],z[n+1:]))
        #compute conditional probabilities of z
        z_conditionals, zvals = conditprob_zj(zs_minusn,x_minusn,x_n,hyperbeta,hypergamma,hyperalpha)
        #workaround for now
        z[n] = sample_z_log(zvals, z_conditionals)
#update parameter values
p = update_bernoullip(z,data)

save_obj(p, "SRR1947693_bernoulliparameters.dict")

