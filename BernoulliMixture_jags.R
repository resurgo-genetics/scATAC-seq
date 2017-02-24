library(rjags)
library(RColorBrewer)
data = read.delim("/Users/cass/Documents/AppliedMachineLearning/HW1/train.csv",header=T,sep=",",nrow=1000)
labels = data[,1]
data[data>=1] = 1
data = sapply(data,as.numeric)

model_string2 <- "
var z[1000], p[10,784], pi[10]
model{
for(i in 1:1000){
z[i] ~ dcat(pi)
for(j in 1:784){
y[i,j] ~ dbern(p[z[i],j])}}

#priors
pi ~ ddirch(dp1)
for(i in 1:10){
for(j in 1:784){
p[i,j] ~ dbeta(1,1)
}}
}
"

model2 <- jags.model(textConnection(model_string2),data =list(y=data,dp1=structure(rep(10.0,10), .Dim = c(10, 1))))
update(model2, 100); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model2, variable.names=c("p", "pi","z"), n.iter=100)
mcmc_summary = summary(mcmc_samples)
test = mcmc_summary$statistics
testnames = row.names(test)
p = test[which(startsWith(testnames,"p[9,"))]
write.table(p,"test.out",quote=F,row.names=F,col.names=F,sep="\t")
