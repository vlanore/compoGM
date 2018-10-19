# Rscript effectiveSize.R <input> > <output>
library(coda)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
d <- read.table(input, header = T, sep = '\t')
r <- raftery.diag(mcmc(d))
write.table(r$resmatrix[,1], file = output, row.names = T, col.names = F, sep = '\t', quote=F)

