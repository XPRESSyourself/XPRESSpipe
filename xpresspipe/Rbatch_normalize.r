#!/usr/bin/env Rscript

# Control batch effects for prep, chips, etc

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("BH", quietly = TRUE)) {install.packages("BH", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("cpp11", quietly = TRUE)) {install.packages("cpp11", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("plogr", quietly = TRUE)) {install.packages("plogr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("XML", quietly = TRUE)) {install.packages("XML", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("dbplyr", quietly = TRUE)) {install.packages("dbplyr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("data.table", quietly = TRUE)) {install.packages("data.table", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("pamr", quietly = TRUE)) {install.packages("pamr", repos = "http://cran.us.r-project.org", quiet = TRUE)}

if ("sva" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("sva", dependencies=c("Depends", "Imports"), quiet = TRUE)
} else {
  print("sva package already installed")
}
if ("bladderbatch" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("bladderbatch", dependencies=c("Depends", "Imports"), quiet = TRUE)
} else {
  print("bladderbatch package already installed")
}
if ("limma" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("limma", dependencies=c("Depends", "Imports"), quiet = TRUE)
} else {
  print("limma package already installed")
}

library(sva)
library(bladderbatch)
library(pamr)
library(limma)

# Get arguments
# args[1] = dataframe
# args[2] = batch info
# args[5] = output file name
args = commandArgs(trailingOnly=TRUE)

# Import expression matrix and sample info
data = read.csv(args[1], sep='\t', row.names=1)
data_m = as.matrix(data)
info = read.csv(args[2], sep='\t')

# Identify batch effect to normalize on
names(info) <- tolower(names(info))
batch = info$batch

# Prepare model
modcombat = model.matrix(~1, data=info)

# Perform combat normalization on known batch parameter
combat_data = ComBat(dat=data_m, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.table(as.data.frame(combat_data), file=args[3], sep='\t', col.names=T, row.names=T)