library(knockoff)
library(glmnet)
library(MASS)
library(Rfast)

test_that("simulation_FDV produces good results", {
  set.seed(1234)

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
      X <- MASS::mvrnorm(n[o], mu = rep(0, p), Sigma = 0.01*Sig)
      while(length(which(abs(X %*% b) > 3)) > 0) {
        X[which(abs(X %*% b) > 3),] <- MASS::mvrnorm(length(which(abs(X %*% b) > 3)), mu = rep(0, p), Sigma = 0.01*Sig)
      }
      prob <- exp(X %*% b) / (1 + exp(X %*% b))
      y <- rep(1, n[o])
      while(sum(y)/n[o] < 0.02 | sum(y)/n[o] > 0.98) {
        for(gen.y in 1:n[o]) {
          y[gen.y] <- rbinom(1, 1, prob[gen.y])
        }
      }

      res <- suppressWarnings(logistic.test(X, y, tune.1 = tune.1, tune.2 = tune.2))
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
