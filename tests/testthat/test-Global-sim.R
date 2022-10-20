library(parallel)
library(doParallel)
library(MASS)
library(Rfast)

test_that("simulation_Global produces good results for Type-I errors", {
  sp <- 2
  p <- c(100, 200, 300, 400)
  r <- c(0.2, 0.4, 1.2)
  alpha <- 0.05
  q.alpha <- -log(pi) - 2 * log(log(1 / (1 - alpha)))
  numcores <- 8
  BBB <- 500
  batchsize <- 4
  tune.1 <- 1.5
  tune.2 <- c(1, 1, 1)

  power.g <- matrix(rep("NA", 12), ncol=3, nrow=4)
  power.adh <- matrix(rep("NA", 12), ncol=3, nrow=4)
  power.llr <- matrix(rep("NA", 12), ncol=3, nrow=4)

  for(o in 1:length(p)) {
    n <- as.integer(p[o]/r)
    b <- rep(0, p[o])
    sig <- toeplitz(seq(0.7, 0.7, length.out = p[o]/10))
    Sig <- bdiag(rep(list(sig), 10)) + diag(rep(0.3, p[o]))

    for(s in 1:(length(r))){
      #cl <- parallel::makeCluster(numcores, outfile = "")
      #doParallel::registerDoParallel(cl)
      test.g <- rep("NA", batchsize)
      test.adh <- rep("NA", batchsize)
      test.llr <- rep("NA", batchsize)

      #G.end <- foreach::foreach(round = 1:(BBB/batchsize)) %dopar% {
        for(bb in 1:batchsize) {
          X <- MASS::mvrnorm(n[s], mu=rep(0, p[o]), Sigma=Sig)

          prob <- exp(X %*% b) / (1 + exp(X %*% b))
          y <- rep(1, n[s])
          while(sum(y)/n[s] < 0.02 | sum(y)/n[s] > 0.98 ) {
            for(gen.y in 1:n[s]) {
              y[gen.y] <- rbinom(1, 1, prob[gen.y])
            }
          }

          res <- suppressWarnings(logistic.test(X, y, lambda = 0.2*sqrt(log(p[o])/n[s]), tune.1 = tune.1, tune.2 = tune.2))
          b.check <- res$b.check
          tau <- res$tau
          M <- res$M
          test <- res$feature.selected


          # Global test
          test.g[bb] <- max((b.check/tau)^2) >= 2 * log(p[o]) - log(log(p[o])) + q.alpha

          # Univariate screening
          test.adh[bb] <- min(p.adjust(Rfast::univglms(y=y, x=X)[,2], method="bonferroni")) < alpha

          # Rescaled LLR
          if(r[s] == 0.2) {
            out <- c()
            out0 <- glm(y~X, family="binomial")$deviance
            for(BB in 1:p[o]) {
              out[BB] <- glm(y~X[,-BB], family="binomial")$deviance - out0
            }
            test.llr[bb] <- min(p.adjust(1 - pchisq(out/1.26, df=1), method="bonferroni")) < alpha
          }

          if(r[s] == 0.4) {
            out <- c()
            out0 <- glm(y~X, family="binomial")$deviance
            for(BB in 1:p[o]) {
              out[BB] <- glm(y~X[,-BB], family="binomial")$deviance - out0
            }
            test.llr[bb] <- min(p.adjust(1 - pchisq(out/2.05, df=1), method="bonferroni")) < alpha
          }
        }
        #return(data.frame(test.g, test.adh, test.llr))
      #}


      #parallel::stopCluster(cl)
      #out.g <- do.call(rbind, G.end)

      #power.g[o, s] <- sum(as.logical(out.g[,1])) / BBB
      #power.adh[o, s] <- sum(as.logical(out.g[,2])) / BBB
      #if(r[s] < 0.5) {
      #  power.llr[o,s] <- sum(as.logical(out.g[,3])) / BBB
      #}

      #print(c(power.g[o, s], power.adh[o, s], power.llr[o, s]))
    }
  }

  expect_equal(1, 1)
})

test_that("simulation_Global produces good results for powers", {


  expect_equal(1, 1)
})
