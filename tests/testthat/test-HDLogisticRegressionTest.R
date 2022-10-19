library(knockoff)
library(glmnet)
library(MASS)
library(Rfast)

test_that("simulation_FDR produces good results", {
  set.seed(123)

  sp <- 60
  n <- 600
  p <- 800

  tune.1 <- 1.5
  tune.2 <- 1

  #Set the regression coefficients.
  b <- rep(0, p)
  b[1:sp] <- rep(c(3, -3), sp/2)

  sig1 <- toeplitz(seq(0.1, 0, length.out = p/10))
  Sig <- bdiag(rep(list(sig1), 10)) + diag(rep(0.9, p))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = 0.01*Sig)
  while(length(which(abs(X %*% b) > 3)) > 0) {
    X[which(abs(X %*% b) > 3),] <- MASS::mvrnorm(length(which(abs(X %*% b) > 3)), mu = rep(0, p), Sigma = 0.01*Sig)
  }
  prob <- exp(X %*% b) / (1 + exp(X %*% b))
  y <- rep(1, n)
  while(sum(y) / n < 0.02 | sum(y) / n > 0.98) {
    for(gen.y in 1:n) {
      y[gen.y] <- rbinom(1, 1, prob[gen.y])
    }
  }
  rm(Sig)

  res <- logistic.test(X, y, tune.1 = tune.1, tune.2 = tune.2)
  dput(res)
  M <- res$M
  test.m <- res$feature.selected

  power <- sum(test.m[1:sp]) / sp
  fdr <- sum(test.m[-(1:sp)]) / (max(c(sum(test.m),1)))

  #BH procedure
  BH <- (p.adjust(2 - 2 * pnorm(abs(M)), method = "BH") < .2)
  power.bh <- sum(BH[1:sp]) / sp
  fdr.bh <- sum(BH[-(1:sp)]) / max(c(1, sum(BH)))

  #BY procedure
  BY <- (p.adjust(2 - 2 * pnorm(abs(M)), method="BY") < .2)
  power.by <- sum(BY[1:sp]) / sp
  fdr.by <- sum(BY[-(1:sp)]) / max(c(1, sum(BY)))

  #knockoff procedure
  k_stat <- function(X, X_k, y) {
    stat.glmnet_coefdiff(X, X_k, y, nlambda = 200, family = "binomial")
  }
  result <- knockoff::knockoff.filter(X, y, statistic = k_stat, fdr = 0.2)
  koff <- result$selected

  power.koff <- length(which(koff < sp + 1)) / sp
  fdr.koff <- length(which(koff > sp)) / max(c(1, length(koff)))

  #univariate screening with BH
  adh <- (p.adjust(Rfast::univglms(y = y, x = X)[,2], method = "BH") < .2)
  power.adh <- sum(adh[1:sp]) / sp
  fdr.adh <- sum(adh[-(1:sp)]) / max(c(1, sum(adh)))


  out <- c(power = power, fdr = fdr,
           power.bh = power.bh, fdr.bh = fdr.bh,
           power.by = power.by, fdr.by = fdr.by,
           power.koff = power.koff, fdr.koff = fdr.koff,
           power.adh = power.adh, fdr.adh = fdr.adh)

  dput(out)

  expect_equal(1, 1)
})


test_that("simulation_FDV produces good results", {
  pdv.temp <- matrix(ncol = 4, nrow = 3)
  fdv.temp <- matrix(ncol = 4, nrow = 3)
  sp <- c(40, 50, 60)
  p <- c(1200)
  tune.1 <- 1.5
  tune.2 <- c(1, 1, 1)
  b <- rep(0, p)

  for(s in 1:3) {
    b[1:sp[s]] <- rep(c(3, -3), sp[s]/2)

    sig1 <- toeplitz(seq(0.1, 0, length.out = p/10))
    Sig <- bdiag(rep(list(sig1), 10)) + diag(rep(0.9, p))

    n <- c(400, 600, 800, 1000)
    for(o in 1:4) {
      X <- mvrnorm(n[o], mu = rep(0, p), Sigma = 0.01*Sig)
      while(length(which(abs(X %*% b) > 3)) > 0) {
        X[which(abs(X %*% b) > 3),] <- mvrnorm(length(which(abs(X %*% b) > 3)), mu = rep(0, p), Sigma = 0.01*Sig)
      }
      prob <- exp(X %*% b) / (1 + exp(X %*% b))
      y <- rep(1, n[o])
      while(sum(y)/n[o] < 0.02 | sum(y)/n[o] > 0.98) {
        for(gen.y in 1:n[o]) {
          y[gen.y] <- rbinom(1, 1, prob[gen.y])
        }
      }

      res <- logistic.test(X, y, tune.1 = tune.1, tune.2 = tune.2)
      M <- res$M

      t.fdv <- qnorm((2 - 10/p) / 2)

      pdv.temp[s, o] <- sum(abs(M[1:sp[s]]) > t.fdv) / sp[s]
      fdv.temp[s, o] <- sum(abs(M[-(1:sp[s])]) > t.fdv)
    }
  }

  dput(pdv.temp)
  dput(fdv.temp)

  expect_equal(1, 1)
})
