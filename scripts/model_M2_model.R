
model <- "model {
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      K[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- sf[j] * 10 ** log10_q[i,cond[j]] * tau[i,j]
      tau[i,j] ~ dgamma(1/alpha[i], 1/alpha[i]) T(0.00001,)
    }
  }

  for(i in 1:ngenes) {
    for(k in 1:ncond) {
      log10_q[i, k] ~ dnorm(3, 1 / (1.5 ** 2))
    }
  }

  for(i in 1:ngenes) {
    log10_alpha[i] ~ dnorm(-2, 1/2*2)
    alpha[i] <- 10 ** log10_alpha[i]
  }
}
"
