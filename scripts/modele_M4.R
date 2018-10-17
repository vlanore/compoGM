## run with Rscript model_M2.R <input-directory> <output-directory>

model <- "model {
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      K[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- sf[j] * q[i,cond[j]] * tau[i,j]
      tau[i,j] ~ dgamma(theta[i,cond[j]], theta[i,cond[j]])
    }
  }

  for(i in 1:ngenes) {
    for(k in 1:ncond) {
      theta[i, k] <- thetamod[i] * (1 / 10 ** (a0 + a1 / q[i, k]))
      q[i, k] <- 10 ** log10_q[i, k]
      log10_q[i, k] ~ dnorm(2, 1 / (2 ** 2))
    }
  }

  for(i in 1:ngenes) {
    thetamod[i] ~ dgamma(1,1)
  }
  a0 ~ dlnorm(-2, 1/2**2)
  a1 ~ dlnorm(0, 1/2**2)
}
"
