#!/usr/bin/env Rscript

# Measure library complexity of RNA-seq sample

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org")}

if ("Rsubread" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("Rsubread", dependencies=TRUE)
} else {
  print("Rsubread package already installed")
}
if ("dupRadar" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("dupRadar", dependencies=TRUE)
} else {
  print("dupRadar package already installed")
}

library(dupRadar)

# Get arguments
# args[1] = de-duplicated BAM file
# args[2] = GTF reference file (full)
# args[3] = is library paired-end? (TRUE or FALSE)
# args[4] = number of threads
# args[5] = output file for metrics
args = commandArgs(trailingOnly=TRUE)

# Set parameters
BAM_NO_DUPS <- args[1]
GTF <- args[2]
STRANDED <- 0 # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
PAIRED <- args[3] # True if paired, or else will assume false
THREADS <- as.numeric(args[4])
OUTPUT <- args[5] # Path and filename with .txt extension

if (PAIRED == 'True') {
  PAIRED = TRUE} else {
  PAIRED = FALSE}

# Duplication rate analysis
dm <- analyzeDuprates(BAM_NO_DUPS, GTF, STRANDED, PAIRED, THREADS)

## Save metrics
write.table(as.data.frame(dm), file=OUTPUT, sep='\t', col.names=NA)
