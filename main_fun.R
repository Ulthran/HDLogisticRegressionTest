#### Input:
#x: n x p matrix of covariates
#y: n dimensional vector.
#nfolds: number of folds used for cross-validation of Lasso/logistic Lasso, default = 5
#lambda: tuning parameter for the logistic lasso, default = NA
#tune.1, tune.2 are tuning parameters in the nodewise lasso. Should be values between 0 and 1
#fdr: false discovery rate for multiple testing


logistic.test <- function(x,y, nfolds=5, lambda = 0, 
                          tune.1 = 1.9, tune.2 = 1.01, 
                          intercept = F, fdr = 0.05){
  f = function(x){
    exp(x)/(1+exp(x))
  }
  p = dim(x)[2]
  n = dim(x)[1]
  if(lambda == 0){
    logistic.cv = cv.glmnet(x = x, y = y, family = "binomial", alpha = 1,  
                            intercept=intercept, nfolds = nfolds, type.measure = "class")
    lambda = logistic.cv$lambda.min
  }
  my.logistic.fit = glmnet(x = x, y = y, family = "binomial", alpha = 1,  
                           intercept=F, lambda = lambda)
  b.hat = coef(my.logistic.fit)
  print("Estimation Succeeded!")
  
  W.n1 = c(exp(x%*% as.matrix(b.hat)[-1,])/(1+exp(x%*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try = matrix(ncol = p,nrow =5)
  tau.try = matrix(ncol = p,nrow = 5)
  
  V = matrix(ncol=n, nrow = p)
  tau = c()
  for(i in 1:p){
    nodewise.try = glmnet(x= x[,-i], y = x[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for(lambda.i in 1:5){
      V[i,] = x[,i]-x[,-i]%*%nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%x[,-i])/sqrt(sum((V[i,])^2*W.n1))))
      tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% x[,i])
    }
    zeta0 = sqrt(2*log(p))
    if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
    lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
    tau[i] = tau.try[lambda.chosen,i]
    lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
    tau[i] = tau.try[lambda.chosen,i]
    V[i,] = x[,i]-x[,-i]%*%nodewise.try$beta[,lambda.chosen]
  }
  
  
  V2 = t((t(V)*W.n1))
  #debaised estimator
  b.check = c()
  for(j in 1:p){
    b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-f(x %*% as.matrix(b.hat)[-1])))/(V[j,] %*% x[,j])
  }
  print("Bias Corrected!")
  
  
  M = b.check/tau
  
  cutoff = function(x) p*(1-1*pchisq(x^2, df=1))/max(1,sum(abs(M) > x)) - fdr
  t = sqrt(2*log(p))
  gd= seq(0, sqrt(2*log(p)-2*log(log(p))), 0.001)
  for(k in 1:length(gd)){
    if(cutoff(gd[k])<0) t = min(gd[k],t)
  }
  test = c()
  test[abs(M)>t] = 1
  test[abs(M)<t] = 0
  
  return(list(b.check = b.check, M = M, 
              g.p.value = 1-exp(-1/sqrt(pi)*exp(-max((b.check/tau)^2)/2)),
              feature.selected = test))
}


#Output

#b.check: vector of debiased regression coefficients.
#M: vector of standardized statistics.
#g.p.value: p-value for the global testing.
#feature.selected: indicator of selected features.

