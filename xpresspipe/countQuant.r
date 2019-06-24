#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org")}

if ("tximport" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("tximport", dependencies=TRUE)
} else {
  print("tximport package already installed")
}

library(tximport)

# Get arguments
# args[1] = dataframe_path
# args[2] = quantification_type
# args[3] = output file name
args = commandArgs(trailingOnly=TRUE)

# Import counts_data
quant_table <- read.table(args[1],sep='\t',header=TRUE,row.names=1)

txi <- tximport(quant_table, type = args[2], tx2gene = tx2gene)

# Write output to new file
write.table(as.data.frame(txi$counts), file=args[3], sep='\t', col.names=T, row.names=T)
