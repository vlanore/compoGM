library(rjags)

source('modele_M3.R')
args <- commandArgs(trailingOnly=TRUE)
dest  <- args[1]

counts_file <- paste0(dest, "/counts.tsv")
sample_file <- paste0(dest, "/samples.tsv")
size_factor_file <- paste0(dest, "/size_factors.tsv")

ngenes <- 100
nsamples <- 2 * 10
ncond <- 2
cond <- c(rep(1, nsamples / 2),
          rep(2, nsamples / 2))
sf <- rgamma(nsamples, 100, 100)

data4jags <- function() {
    list(
      sf = sf,
      ngenes = ngenes,
      nsamples = nsamples,
      ncond = ncond,
      cond = cond
    )
}

run_mcmc <- function() {
    mcmc <- jags.model(textConnection(model),
                       data = data4jags(),
                       n.chains = 1)
    sample <- coda.samples(mcmc,
                           variable.names = c("K"),
                           n.iter = 1)
    sample
}

save_sim <- function(sample) {
  m <- matrix(data=(sample[[1]]), nrow = ngenes, ncol = nsamples, byrow=F)
  gene_names <- data.frame(paste0("g", 1:ngenes))
  sample_names <- paste0("s", 1:nsamples)
  d <- cbind(gene_names, data.frame(m))
  colnames(d) <- c("gene", sample_names)
  system(paste0("mkdir -p ", dest))
  write.table(sf, file = size_factor_file, sep = "\t", quote = F, row.names = F)
  write.table(d, file = counts_file, sep = "\t", quote = F, row.names = F)
  write.table(data.frame(sample = sample_names, condition = cond), file = sample_file, sep = "\t", quote = F, row.names = F)
}

trace <- run_mcmc()
save_sim(trace)
