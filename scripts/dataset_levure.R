library(DESeq2)

directory <- "."
sampleFiles <- grep("gbgout",list.files(directory),value=TRUE)
sampleCondition <- sub("([^_]+).*","\\1",sampleFiles)
sampleName <- sub("([^_]+_[^_]+).*","\\1",sampleFiles)

kept <- !(sampleName %in% c("Snf2_rep06","Snf2_rep13","Snf2_rep35",
                            "WT_rep34","WT_rep36","WT_rep21",
                            "WT_rep22","WT_rep25","WT_rep28"))

counts <- lapply(sampleFiles[kept], function(fn) {
  head(read.table(fn)[[2]], n=-5)
})
counts <- data.frame(counts)
gene_names <- head(read.table(sampleFiles[1])[[1]], n=-5)
colnames(counts) <- sampleName[kept]
count_table <- cbind(gene = gene_names, counts)

sample_table <- data.frame(sample = sampleName[kept],
                           condition = sampleCondition[kept])

sampleTable <- data.frame(sampleName = sampleName,
                          fileName = paste(directory, sampleFiles, sep="/"),
                          condition = sampleCondition)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable[kept,],
                                  directory = ".",
                                  design = ~ condition)
dds <- estimateSizeFactors(dds)

size_factors <- data.frame(sample = sampleName[kept], size_factor = sizeFactors(dds))

write.table(file = "counts.tsv", count_table, quote = F, row.names = F, sep = "\t")
write.table(file = "samples.tsv", sample_table, quote = F, row.names = F, sep = "\t")
write.table(file = "size_factors.tsv", size_factors, quote = F, row.names = F, sep = "\t")

