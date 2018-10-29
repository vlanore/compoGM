
model <- "model {
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      K[i,j] ~ dpois(10 ** log10_lambda[i, cond[j]])
    }
  }

  for(i in 1:ngenes) {
    for(k in 1:ncond) {
      log10_lambda[i, k] ~ dnorm(3, 1 / (1.5 ** 2))
    }
  }

}
"
