import pymc3 as pm
import numpy as np

train_csv = open("/Users/cass/Documents/AppliedMachineLearning/HW1/train.csv","r")
lines = train_csv.readlines()[1:]
data = np.zeros([len(lines),len(lines[0].split(','))-1])

truth = []; i=0
for line in lines:
    data[i] = map(float,line.split(',')[1:])
    truth.append(line[0])
    i+=1
truth = map(int,truth)

#just take a subset
data = data[0:1000]
truth = truth[0:1000]

data[np.where(data>=1)]=1


np.random.seed(12345)
alphaprime=10
nclusters = 10
ncells = data.shape[0]
nsites = data.shape[1]

#without scaling
model = pm.Model()
with model:
    pi = pm.Dirichlet('pi', a=np.array([alphaprime]*nclusters),shape=nclusters)
    # Define priors
    pk = pm.Beta('pk', 1,1,shape=(nclusters,nsites))
    z = pm.Categorical("z",p=pi,shape=ncells)
    # Define likelihood
    likelihood = pm.Bernoulli('likelihood',p=pk[z],observed=data,shape=(ncells))


with model:
    step1 = pm.Metropolis(vars=[pk, pi, alpha])
    step2 = pm.ElemwiseCategorical(vars=[category], values=[0, 1, 2])
    tr = pm.sample(100, step=[step1, step2])

traceplot(trace)





#attempt 1: specified as in BISCUIT
model = pm.Model()
with model: # model specifications in PyMC3 are wrapped in a with-statement
    pi = pm.Dirichlet('p', a=np.array([alphaprime]*nclusters))
    # Define priors
    pk = Beta('pk', 1,1,shape=nclusters)
    alpha = Beta('alpha',1,.1,shape=ncells)
    z = Categorical("z",p=pi,shape=ncells)
    # Define likelihood
    likelihood =Bernoulli('y', p=pk[z]*alpha[y],observed=data)
    step1 = pm.Metropolis(vars=[pk, pi, alpha])
    step2 = pm.ElemwiseCategorical(vars=[category], values=[0, 1, 2])
    tr = pm.sample(100, step=[step1, step2])

traceplot(trace)



#attempt 2: specified in PYMC3 examples
def stick_breaking(beta):
    portion_remaining = tt.concatenate([[1], tt.extra_ops.cumprod(1 - beta)[:-1]])
    return beta * portion_remaining

with pm.Model() as model:
    alphaprime = pm.Gamma('alpha', 1., 1.)
    beta = pm.Beta('beta', 1., alphaprime, shape=K)
    w = pm.Deterministic('w', stick_breaking(beta))
    pk = Beta('pk', 1,1,shape=nclusters)
    alpha = Beta('alphaprime',1,.1,shape=ncells)
    likelihood =Bernoulli('y', p=pk[w]*alpha[y],observed=y)

                           