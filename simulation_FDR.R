#The simulations were ran in a parallel computing platform. The following codes consist of single round of simulation with specificed n, s, and p.

args= commandArgs(trailingOnly=TRUE)

library(knockoff)
library(glmnet)
library(doMC)
library(Rfast)
library(MASS)
library(purrr)


sp=60
n=600
p=800

tune.1=1.5
tune.2=1
f = function(x){
  exp(x)/(1+exp(x))
}


#Set the regression coefficients.
b = rep(0,p)
b[1:sp] = rep(c(3,-3), sp/2)


power = 0
fdr = 0
power.bh = 0
fdr.bh = 0
power.by = 0
fdr.by = 0
power.koff=0
fdr.koff=0
power.adh=0
fdr.adh=0
test.m = c()
sig1 = toeplitz(seq(0.1,0,length.out = p/10))
Sig = bdiag(rep(list(sig1),10))+diag(rep(0.9,p))
#Sig=diag(rep(1,p))
X=mvrnorm(n, mu=rep(0,p), Sigma=0.01*Sig)
while(length(which(abs(X%*% b)>3))>0){
  X[which(abs(X%*%b)>3),]=mvrnorm(length(which(abs(X%*% b)>3)), mu=rep(0,p), Sigma=0.01*Sig)
}
prob = exp(X %*% b)/(1+exp(X %*% b))
y = rep(1,n)
while(sum(y)/n<0.02 | sum(y)/n>0.98 ){
  for(gen.y in 1:n){
    y[gen.y]=rbinom(1,1,prob[gen.y])
  }
}
rm(Sig)

####### lasso estimation

my.logistic.fit = glmnet(x = X, y = y, family = "binomial", alpha = 1,  intercept=F,
                         lambda = 0.01*sqrt(log(p)/n), standardize=F)
b.hat = coef(my.logistic.fit)

####### score vector

W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
zeta.try = matrix(ncol = p,nrow =5)
tau.try = matrix(ncol = p,nrow = 5)

V = matrix(ncol=n, nrow = p)
tau = c()
for(i in 1:p){
  nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
  for(lambda.i in 1:5){
    V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
    zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(sum((V[i,])^2*W.n1))))
    tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% X[,i])
  }
  zeta0 = sqrt(2*log(p))
  if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
  lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
  tau[i] = tau.try[lambda.chosen,i]
  lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
  tau[i] = tau.try[lambda.chosen,i]
  V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
}


V2 = t((t(V)*W.n1))
#debaised estimator
b.check = c()
for(j in 1:p){
  b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
}

#### multiple test
M = b.check/tau

cutoff = function(x, alpha) p*(1-1*pchisq(x^2, df=1))/max(1,sum(abs(M) > x)) - alpha

t = sqrt(2*log(p))
x= seq(0, sqrt(2*log(p)-2*log(log(p))), 0.001)
for(k in 1:length(x)){
  if(cutoff(x[k], 0.2)<0) t = min(x[k],t)
}


test.m[abs(M)>t] = 1
test.m[abs(M)<t] = 0

power = sum(test.m[1:sp])/sp
fdr = sum(test.m[-(1:sp)])/(max(c(sum(test.m),1)))

#BH procedure
BH = (p.adjust(2-2*pnorm(abs(M)),method="BH")<.2)
power.bh =sum(BH[1:sp])/sp
fdr.bh = sum(BH[-(1:sp)])/max(c(1,sum(BH)))

#BY procedure
BY = (p.adjust(2-2*pnorm(abs(M)),method="BY")<.2)
power.by = sum(BY[1:sp])/sp
fdr.by = sum(BY[-(1:sp)])/max(c(1,sum(BY)))  

#knockoff procedure
k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "binomial")
result = knockoff.filter(X, y, statistic=k_stat, fdr=0.2)
koff=result$selected

power.koff = length(which(koff<sp+1))/sp
fdr.koff = length(which(koff >sp))/max(c(1,length(koff)))  

#univariate screening with BH
adh = (p.adjust(univglms(y=y,x=X)[,2],method="BH")<.2)
power.adh = sum(adh[1:sp])/sp
fdr.adh = sum(adh[-(1:sp)])/max(c(1,sum(adh)))


out=c(power, fdr, power.bh, fdr.bh, power.by, fdr.by, power.koff, fdr.koff, power.adh, fdr.adh)

save(out, file= args[1])