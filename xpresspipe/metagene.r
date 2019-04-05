#!/usr/bin/env Rscript

# Measure metagene coverage of RNA-seq sample

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("metagene", version = "3.8")
library(metagene)
