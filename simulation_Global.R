#The following codes consist of simulations of global testing with sparsity sp=2.
library(MASS)
library(purrr)
library(glmnet)
library(Rfast)
library(parallel)
library(doParallel)
library(foreach)


#####Global Testing, Type I Errors

sp=2

p=c(100,200,300,400)
r = c(0.2,0.4,1.2)
alpha= 0.05
q.alpha = -log(pi)-2*log(log(1/(1-alpha)))
numcores=8
BBB=500
batchsize=4

tune.1=1.5
tune.2=c(1,1,1)
f = function(x){
  exp(x)/(1+exp(x))
}


power.g = matrix(rep("NA",12),ncol=3, nrow=4)
power.adh = matrix(rep("NA",12),ncol=3, nrow=4)
power.llr= matrix(rep("NA",12),ncol=3, nrow=4)

for(o in 1:length(p)){
  n=as.integer(p[o]/r)
  b = rep(0,p[o])
  sig = toeplitz(seq(0.7,0.7,length.out = p[o]/10))
  Sig = bdiag(rep(list(sig),10))+diag(rep(0.3,p[o])) 
  
  for(s in 1:(length(r))){
    
    cl <- parallel::makeCluster(numcores, outfile = "")
    registerDoParallel(cl)
    test.g = rep("NA",batchsize)
    test.adh=rep("NA",batchsize)
    test.llr=rep("NA",batchsize)
    
    G.end=foreach(round=1:(BBB/batchsize)) %dopar% {
      library(MASS)
      library(purrr)
      library(glmnet)
      library(Rfast)
      for(bb in 1:batchsize){
        X = mvrnorm(n[s], mu=rep(0,p[o]), Sigma=Sig)
        
        prob = exp(X %*% b)/(1+exp(X %*% b))
        y = rep(1,n[s])
        while(sum(y)/n[s]<0.02 | sum(y)/n[s]>0.98 ){
          for(gen.y in 1:n[s]){
            y[gen.y]=rbinom(1,1,prob[gen.y])
          }
        }
        
        ####### biased estimation
        
        my.logistic.fit = glmnet(x = X, y = y, family = "binomial", alpha = 1,  intercept=F, lambda = 0.2*sqrt(log(p[o])/n[s]))
        b.hat = coef(my.logistic.fit)
        
        ####### score vector
        W.n1 = diag(c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1))
        zeta.try = matrix(ncol = p[o],nrow =5)
        tau.try = matrix(ncol = p[o],nrow = 5)
        
        V = matrix(ncol=n[s], nrow = p[o])
        tau = c()
        for(i in 1:p[o]){
          nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5)
          for(lambda.i in 1:5){
            V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
            zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(V[i,]%*%W.n1%*%V[i,])))
            tau.try[lambda.i,i] = sqrt(V[i,]%*%W.n1%*%V[i,])/(V[i,]%*% X[,i])
          }
          zeta0 = sqrt(2*log(p[o]))
          if(min(zeta.try[,i])>sqrt(2*log(p[o]))) zeta0 = tune.1*min(zeta.try[,i])
          lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
          tau[i] = tau.try[lambda.chosen,i]
          lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2[s]*tau[i]],decreasing = F)[1]
          tau[i] = tau.try[lambda.chosen,i]
          V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
        }
        
        
        
        V2 = V%*%W.n1
        #debaised estimator
        b.check = c()
        for(j in 1:p[o]){
          b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
        }
        
        
        #### global test
        test.g[bb] = max((b.check/tau)^2) >= 2*log(p[o])-log(log(p[o]))+q.alpha
        
        
        #univariate screening
        test.adh[bb] = min(p.adjust(univglms(y=y,x=X)[,2],method="bonferroni"))<alpha
        
        #rescaled LLR
        if(r[s]==0.2){
          out=c()
          out0 = glm(y~X, family="binomial")$deviance
          for(BB in 1:p[o]){
            out[BB]=glm(y~X[,-BB], family="binomial")$deviance-out0
          }
          test.llr[bb] = min(p.adjust(1-pchisq(out/1.26,df=1),method="bonferroni"))<alpha
        }
        
        if(r[s]==0.4){
          out=c()
          out0 = glm(y~X, family="binomial")$deviance
          for(BB in 1:p[o]){
            out[BB]=glm(y~X[,-BB], family="binomial")$deviance-out0
          }
          test.llr[bb] = min(p.adjust(1-pchisq(out/2.05,df=1),method="bonferroni"))<alpha
        }
        
        
      }
      return(data.frame(test.g,test.adh,test.llr)) 
    }
    
    
    stopCluster(cl)
    out.g= do.call(rbind, G.end)
    
    power.g[o,s] = sum(as.logical(out.g[,1]))/BBB
    power.adh[o,s]=sum(as.logical(out.g[,2]))/BBB
    if(r[s]<0.5) power.llr[o,s] =sum(as.logical(out.g[,3]))/BBB
    
    print(c(power.g[o,s],power.adh[o,s],power.llr[o,s]))
  }
  
}

save(power.g,power.adh,power.llr, file="Global_Err_sp2.RData")



###### Powers

#####Global Testing

sp=2

p=c(100,200,300,400)
r = c(0.2,0.4,1.2)
alpha= 0.05
q.alpha = -log(pi)-2*log(log(1/(1-alpha)))
numcores=8
BBB=400
batchsize=2

tune.1=1.5
tune.2=c(1,1,1)
f = function(x){
  exp(x)/(1+exp(x))
}


power.g = matrix(rep("NA",12),ncol=3, nrow=4)
power.adh = matrix(rep("NA",12),ncol=3, nrow=4)
power.llr= matrix(rep("NA",12),ncol=3, nrow=4)

for(o in 1:length(p)){
  n=as.integer(p[o]/r)
  b = rep(0,p[o])
  b[1:sp] = rep(c(-0.75,0.75), sp/2)
  sig = toeplitz(seq(0.7,0.7,length.out = p[o]/10))
  Sig = bdiag(rep(list(sig),10))+diag(rep(0.3,p[o])) 
  
  for(s in 1:(length(r))){
    
    cl <- parallel::makeCluster(numcores, outfile = "")
    registerDoParallel(cl)
    test.g = rep("NA",batchsize)
    test.adh=rep("NA",batchsize)
    test.llr=rep("NA",batchsize)
    
    G.end=foreach(round=1:(BBB/batchsize)) %dopar% {
      library(MASS)
      library(purrr)
      library(glmnet)
      library(Rfast)
      for(bb in 1:batchsize){
        X = mvrnorm(n[s], mu=rep(0,p[o]), Sigma=Sig)
        
        prob = exp(X %*% b)/(1+exp(X %*% b))
        y = rep(1,n[s])
        while(sum(y)/n[s]<0.02 | sum(y)/n[s]>0.98 ){
          for(gen.y in 1:n[s]){
            y[gen.y]=rbinom(1,1,prob[gen.y])
          }
        }
        
        ####### biased estimation
        
        
        my.logistic.fit = glmnet(x = X, y = y, family = "binomial", alpha = 1,  intercept=F, lambda = 0.25*sqrt(log(p[o])/n[s]))
        b.hat = coef(my.logistic.fit)
        
        ####### score vector
        W.n1 = diag(c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1))
        zeta.try = matrix(ncol = p[o],nrow =5)
        tau.try = matrix(ncol = p[o],nrow = 5)
        
        V = matrix(ncol=n[s], nrow = p[o])
        tau = c()
        for(i in 1:p[o]){
          nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5)
          for(lambda.i in 1:5){
            V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
            zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(V[i,]%*%W.n1%*%V[i,])))
            tau.try[lambda.i,i] = sqrt(V[i,]%*%W.n1%*%V[i,])/(V[i,]%*% X[,i])
          }
          zeta0 = sqrt(2*log(p[o]))
          if(min(zeta.try[,i])>sqrt(2*log(p[o]))) zeta0 = tune.1*min(zeta.try[,i])
          lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
          tau[i] = tau.try[lambda.chosen,i]
          lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2[s]*tau[i]],decreasing = F)[1]
          tau[i] = tau.try[lambda.chosen,i]
          V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
        }
        
        
        
        V2 = V%*%W.n1
        #debaised estimator
        b.check = c()
        for(j in 1:p[o]){
          b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
        }
        
        
        #### global test
        test.g[bb] = max((b.check/tau)^2) >= 2*log(p[o])-log(log(p[o]))+q.alpha
        test.adh[bb] = min(p.adjust(regression(X,y)[,2],method="bonferroni"))<alpha
        
        if(r[s]==0.2){
          out=c()
          out0 = glm(y~X, family="binomial")$deviance
          for(BB in 1:p[o]){
            out[BB]=glm(y~X[,-BB], family="binomial")$deviance-out0
          }
          test.llr[bb] = min(p.adjust(1-pchisq(out/1.26,df=1),method="bonferroni"))<alpha
        }
        
        if(r[s]==0.4){
          out=c()
          out0 = glm(y~X, family="binomial")$deviance
          for(BB in 1:p[o]){
            out[BB]=glm(y~X[,-BB], family="binomial")$deviance-out0
          }
          test.llr[bb] = min(p.adjust(1-pchisq(out/2.05,df=1),method="bonferroni"))<alpha
        }
        
        
      }
      return(data.frame(test.g,test.adh,test.llr)) 
    }
    
    
    stopCluster(cl)
    out.g= do.call(rbind, G.end)
    
    power.g[o,s] = sum(as.logical(out.g[,1]))/BBB
    power.adh[o,s]=sum(as.logical(out.g[,2]))/BBB
    if(r[s]<0.5) power.llr[o,s] =sum(as.logical(out.g[,3]))/BBB
    
    print(c(power.g[o,s],power.adh[o,s],power.llr[o,s]))
  }
  
}

save(power.g,power.adh,power.llr, file="Global_Pw_sp2.RData")


