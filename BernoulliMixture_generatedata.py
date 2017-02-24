#code to generate data using generative model for infinite bernoulli mixture
import scipy.stats
import numpy as np
import pymc3 as pm
import math

#number of cell types = K
#number of sites = D
#number of cells = N
#proportion of each cluster in the data is a vector pi

K = 5
D = 500
N=500
pi = [.1,.25,.4,.07,.18]
clusters = [0,1,2,3,4]


#1: set separate hyperparameters for each site
beta1 = {}
gamma = {}
for d in range(D):
    beta1[d] = np.random.exponential(5)
    gamma[d] = beta = np.random.exponential(5)

#sample cluster parameters from prior distributions
p = {}
for k in range(K):
    p[k] = [0]*D
    for d in range(D):
        p[k][d] = np.random.beta(beta1[d],gamma[d])

#2: sample cell-specific scaling factors for technical variation
alpha = {}
for n in range(N):
    alpha[n] = np.random.beta(180,75)

#3: for each cell generate data according to the model
data = np.zeros(shape=(500,500))
for n in range(N):
    #choose a cluster
    k = np.random.choice(clusters,p=pi)
    data[n] = np.random.binomial([1]*D,np.multiply(alpha[n],p[k]))

#now run algorithm


#attempt 1: my implementation without scaling factor
for i in range(50):
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


#attempt 2: specified as in BISCUIT
alphaprime=10
model = pm.Model()
with model: # model specifications in PyMC3 are wrapped in a with-statement
    pi1 = pm.Dirichlet('pi', a=[alphaprime]*k)
    # Define priors
    pk = Beta('pk', 1,1,shape=k)
    alpha1 = Beta('alpha',1,.1,shape=N)
    z = Categorical("z",p=pi,shape=N)
    # Define likelihood
    likelihood =Bernoulli('y', p=pk[z]*alpha1[y],observed=data)
    step1 = pm.Metropolis(vars=[pk, pi1, alpha1])
    step2 = pm.ElemwiseCategorical(vars=[z], values=[0, 1, 2])
    tr = pm.sample(10000, step=[step1, step2])
traceplot(trace)

#attempt 3: specified as in BISCUIT
def stick_breaking(beta):
    portion_remaining = tt.concatenate([[1], tt.extra_ops.cumprod(1 - beta)[:-1]])
    return beta * portion_remaining

with pm.Model() as model:
    alphaprime = pm.Gamma('alpha', 1., 1.)
    beta1 = pm.Beta('beta1', 1., alphaprime, shape=K)
    w = pm.Deterministic('w', stick_breaking(beta))
    pk = Beta('pk', 1,1,shape=K)
    alpha = Beta('alphaprime',1,.1,shape=N)
    likelihood =Bernoulli('y', p=pk[w]*alpha[data],observed=data)
    trace = pm.sample(2000, n_init=50000, random_seed=SEED)
traceplot(trace)
