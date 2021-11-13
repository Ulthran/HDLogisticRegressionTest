#The following codes consist of the simulations of FDV/power for a given p, which is set to be 1200 below.
args= commandArgs(trailingOnly=TRUE)

library(knockoff)
library(glmnet)
library(doMC)
library(Rfast)
library(MASS)
library(purrr)

pdv.temp = matrix(ncol=4,nrow=3)
fdv.temp = matrix(ncol=4,nrow=3)

sp=c(40,50,60)

p=c(1200)

tune.1=1.5
tune.2=c(1,1,1)
f = function(x){
  exp(x)/(1+exp(x))
}

b = rep(0,p)

for(s in 1:3){
  
  b[1:sp[s]] = rep(c(3,-3), sp[s]/2)
  
  
  sig1 = toeplitz(seq(0.1,0,length.out = p/10))
  Sig = bdiag(rep(list(sig1),10))+diag(rep(0.9,p))
  
  n=c(400,600,800,1000)
  for(o in 1:100){
    X=mvrnorm(n[o], mu=rep(0,p), Sigma=0.01*Sig)
    while(length(which(abs(X%*% b)>3))>0){
      X[which(abs(X%*%b)>3),]=mvrnorm(length(which(abs(X%*% b)>3)), mu=rep(0,p), Sigma=0.01*Sig)
    }
    prob = exp(X %*% b)/(1+exp(X %*% b))
    y = rep(1,n[o])
    while(sum(y)/n[o]<0.02 | sum(y)/n[o]>0.98 ){
      for(gen.y in 1:n[o]){
        y[gen.y]=rbinom(1,1,prob[gen.y])
      }
    }
    
    
    ####### biased estimation
    
    my.logistic.fit = glmnet(x = X, y = y, family = "binomial", alpha = 1,  intercept=F,
                             lambda = 0.01*sqrt(log(p)/n[o]), standardize=F)
    b.hat = coef(my.logistic.fit)
    
    ####### score vector
    
    W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
    zeta.try = matrix(ncol = p,nrow =5)
    tau.try = matrix(ncol = p,nrow = 5)
    
    V = matrix(ncol=n[o], nrow = p)
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
    
    t.fdv = qnorm((2-10/p)/2)
    
    pdv.temp[s,o] = sum(abs(M[1:sp[s]])>t.fdv)/sp[s]
    fdv.temp[s,o] = sum(abs(M[-(1:sp[s])])>t.fdv)
    
  }
}

save(pdv.temp,fdv.temp, file= args[1])


