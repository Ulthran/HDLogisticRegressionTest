library('glmnet')

#' DESCRIPTION
#'
#' @param x is a n x p matrix of covariates
#' @param y is a n dimensional vector.
#' @param nfolds is the number of folds used for cross-validation of Lasso/logistic Lasso (Default: 5)
#' @param lambda is the tuning parameter for the logistic lasso (Default: 0)
#' @param tune.1 is a tuning parameter in the nodewise lasso. Should be values between 1 and 2 (Default: 1.9)
#' @param tune.2 is a tuning parameter in the nodewise lasso. Should be values between 1 and 2 (Default: 1.01)
#' @param intercept is a boolean SOMETHING (Default: F)
#' @param fdr is the false discovery rate for multiple testing (Default: 0.05)
#'
#' @return is a named list containing \itemize{
#' \item \code{b.check} is a vector of debiased regression coefficients
#' \item \code{M} is a vector of standardized statistics
#' \item \code{g.p.value} is the p-value for the global testing
#' \item \code{feature.selected} is an indicator of selected features
#' }
#'
#' @importFrom stats coef
#' @importFrom stats pchisq

logistic.test <- function(X, y, nfolds = 5, lambda = 0,
                          tune.1 = 1.9, tune.2 = 1.01,
                          intercept = F, fdr = 0.05){
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (lambda == 0) {
    lambda <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", alpha = 1,
                                     intercept=intercept, nfolds = nfolds, type.measure = "class")$lambda.min
  }

  my.logistic.fit <- glmnet::glmnet(x = X, y = y, family = "binomial", alpha = 1,
                                    intercept=F, lambda = lambda)
  b.hat <- coef(my.logistic.fit)
  print("Estimation Succeeded!")

  W.n1 <- c(exp(X %*% as.matrix(b.hat)[-1,])/(1+exp(X %*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try <- matrix(ncol = p, nrow = 5)
  tau.try <- matrix(ncol = p, nrow = 5)

  V <- matrix(ncol = n, nrow = p)
  tau <- c()
  for (i in 1:p) {
    nodewise.try <- glmnet::glmnet(x = X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for (lambda.i in 1:5) {
      V[i,] <- X[,i] - X[,-i] %*% nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i, i] <- max(abs(as.vector(V[i,] %*% X[,-i]) / sqrt(sum((V[i,])^2 * W.n1))))
      tau.try[lambda.i, i] <- sqrt(sum((V[i,])^2 * W.n1)) / (V[i,] %*% X[,i])
    }
    zeta0 <- sqrt(2 * log(p))
    if (min(zeta.try[,i]) > zeta0) {
      zeta0 <- tune.1 * min(zeta.try[,i])
    }
    lambda.chosen <- order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing = T)[1]
    tau[i] <- tau.try[lambda.chosen, i]
    lambda.chosen <- order(nodewise.try$lambda[tau.try[,i] <= tune.2 * tau[i]], decreasing = F)[1]
    tau[i] <- tau.try[lambda.chosen, i]
    V[i,] <- X[,i] - X[,-i] %*% nodewise.try$beta[,lambda.chosen]
  }

  V2 <- t((t(V) * W.n1))
  # Debiased estimator
  b.check <- c()
  for(j in 1:p){
    b.check[j] <- b.hat[j+1] + (V2[j,] %*% (y - exp(X %*% as.matrix(b.hat)[-1]) / (1 + exp(X %*% as.matrix(b.hat)[-1])))) / (V[j,] %*% X[,j])
  }
  print("Bias Corrected!")


  M <- b.check / tau

  cutoff <- function(x) p*(1-1*pchisq(x^2, df=1))/max(1,sum(abs(M) > x)) - fdr
  t <- sqrt(2 * log(p))
  gd <- seq(0, sqrt(2 * log(p) - 2 * log(log(p))), 0.001)
  for(k in 1:length(gd)){
    if(p * (1 - pchisq(gd[k]^2, df = 1)) / max(1, sum(abs(M) > gd[k])) - fdr < 0){
      t <- min(gd[k], t)
    }
  }

  test <- c()
  test[abs(M)>t] <- 1
  test[abs(M)<t] <- 0

  return(list(b.check = b.check, M = M,
              g.p.value = 1-exp(-1/sqrt(pi)*exp(-max((b.check/tau)^2)/2)),
              feature.selected = test))
}



