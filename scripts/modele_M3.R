## run with Rscript model_M2.R <input-directory> <output-directory>

model <- "model {
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      K[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- sf[j] * 10 ** log10_q[i,cond[j]] * tau[i,j]
      tau[i,j] ~ dgamma(1/alpha[i], 1/alpha[i])
    }
  }

  for(i in 1:ngenes) {
    for(k in 1:ncond) {
      log10_q[i, k] ~ dnorm(2, 1 / (2 ** 2))
    }
  }

  for(i in 1:ngenes) {
    alpha[i] <- 10 ** log10_alpha[i]
    log10_alpha[i] ~ dnorm(log(alpha_bar[i]) / log(10), 1 / sigma_alpha ** 2)
    alpha_bar[i] <- a0 + a1 / q_bar[i]
    q_bar[i] <- sum(10 ** log10_q[i,1:ncond]) / ncond
  }
  sigma_alpha ~ dexp(1)
  a0 ~ dlnorm(-2, 1/2**2)
  a1 ~ dlnorm(0, 1/2**2)
}
"
