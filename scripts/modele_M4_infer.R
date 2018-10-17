library(rjags)

source('modele_M4.R')
args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
dest  <- args[2]
n.chains <- 3

counts <- read.table(paste0(input_dir, "/counts.tsv"), header = T, sep = "\t")
samples <- read.table(paste0(input_dir, "/samples.tsv"), header = T, sep = "\t")
size_factors <- read.table(paste0(input_dir, "/size_factors.tsv"), header = T, sep = "\t")

data4jags <- function(genesel, samplesel) {
    list(
      sf = size_factors[samplesel,2],
      ngenes = length(genesel),
      nsamples = length(samplesel),
      K = counts[genesel, 1 + samplesel],
      ncond = length(unique(samples[samplesel,2])),
      cond = as.integer(samples[samplesel,2])
    )
}

run_mcmc <- function(genesel, samplesel) {
    mcmc <- jags.model(textConnection(model),
                       data = data4jags(genesel, samplesel),
                       n.chains = n.chains)
    update(mcmc, 5000)
    sample <- coda.samples(mcmc,
                           variable.names = c("log10_q", "a0", "a1"),
                           n.iter = 5000)
    sample
}

save_trace <- function(trace) {
  system(paste0("mkdir -p ", dest))
  for(i in 1:n.chains) {
    write.table(trace[[i]],file = paste0(dest,"/trace.",i,".tsv"), sep = "\t", quote = F, row.names = F)
  }
}

trace <- run_mcmc(1:4, 1:dim(samples)[1])
save_trace(trace)
