#!/usr/bin/env Rscript

# Measure library complexity of RNA-seq sample

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dupRadar", version = "3.8")
library(dupRadar)

# Get arguments
# args[1] = de-duplicated BAM file
# args[2] = GTF reference file (full)
# args[3] = is library paired-end? (TRUE or FALSE)
# args[4] = number of threads
# args[5] = output file for metrics

# Set parameters
bamDuprm <- args[1]
gtf <- args[2]
stranded <- 0 # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
paired <- args[3]
threads <- args[4]
output <- args[5]

# Duplication rate analysis
dm <- analyzeDuprates(bamDuprm, gtf, stranded, paired, threads)

## Save metrics
write.table(as.data.frame(dm), file=output, sep='\t', col.names=T, row.names=T)
