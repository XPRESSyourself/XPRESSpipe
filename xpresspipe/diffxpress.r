#!/usr/bin/env Rscript

# Control batch effects for prep, chips, etc
library(DESeq2)

# Get arguments
# args[1] = dataframe
# args[2] = sample info
# args[3] = output file name
# args[4] = DE equation
args = commandArgs(trailingOnly=TRUE)

# Import counts_data
count_table <- read.table(args[1],sep='\t',header=TRUE,row.names=1)

# Create conditions dataframe
sample_table <- read.table(text=readLines(args[2], warn = FALSE), header=TRUE, sep='\t')
names(sample_table) <- tolower(names(sample_table))

# Run DESeq2 analysis on data
dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = sample_table,
  design = as.formula(paste('~', tolower(toString(args[4])))))
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]

# Write output to new file
write.table(as.data.frame(resOrdered), file=args[3], sep='\t', col.names=T, row.names=T)
