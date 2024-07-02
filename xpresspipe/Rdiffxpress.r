#!/usr/bin/env Rscript
license <- function() {
  "
  XPRESSpipe
  An alignment and analysis pipeline for RNAseq data
  alias: xpresspipe

  Copyright (C) 2019  Jordan A. Berg
  jordan <dot> berg <at> biochem <dot> utah <dot> edu

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program.  If not, see <https://www.gnu.org/licenses/>.
  "
  }

# Control batch effects for prep, chips, etc
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("BH", quietly = TRUE)) {install.packages("BH", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("cpp11", quietly = TRUE)) {install.packages("cpp11", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("plogr", quietly = TRUE)) {install.packages("plogr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("XML", quietly = TRUE)) {install.packages("XML", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("dbplyr", quietly = TRUE)) {install.packages("dbplyr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("data.table", quietly = TRUE)) {install.packages("data.table", repos = "http://cran.us.r-project.org", quiet = TRUE)}

if ("DESeq2" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("DESeq2", dependencies=c("Depends", "Imports"), quiet = TRUE)
} else {
  print("DESeq2 package already installed")
}

if ("apeglm" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("apeglm", dependencies=c("Depends", "Imports"), quiet = TRUE)
} else {
  print("apeglm package already installed")
}

library(DESeq2)
library(apeglm)

# Get arguments
# args[1] = dataframe
# args[2] = sample info
# args[3] = output file name
# args[4] = DE equation
args = commandArgs(trailingOnly=TRUE)

# Set parameters
DATAFRAME <- args[1] # Directory containing
INFO <- args[2]
OUTPUT <- args[3] # Path and filename with .txt extension
EQUATION <- args[4]
SHRINKAGE <- args[5]

# Import counts_data
count_table <- read.table(DATAFRAME, sep = '\t', header = TRUE, row.names = 1, check.names=F)

# Create conditions dataframe
sample_table <- read.table(text = readLines(INFO, warn = FALSE), header = TRUE, sep = '\t')
names(sample_table) <- tolower(names(sample_table))

# Convert conditions in design to factor levels
items <- strsplit(gsub("[^[:alnum:] ]", " ", tolower(toString(EQUATION))), " +")[[1]]
uniq_items <- unique(items)

for ( i in uniq_items) {
    sample_table[[i]] <- factor(
      sample_table[[i]], levels = sort(
        unique(
          unlist(
            lapply(
              sample_table[[i]], toString
            )
          )
        )
      )
    )
}

# Run DESeq2 analysis on data
dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = sample_table,
  design = as.formula(paste('~', tolower(toString(EQUATION)))))
dds <- DESeq(dds)

if (SHRINKAGE == TRUE) {
  shrinkage_var <- tail(resultsNames(dds), n=1)
  print(paste('Using factor ', shrinkage_var, ' as log fold change shrinkage coefficient...', sep=''))
  res <- lfcShrink(
    dds,
    coef=toString(shrinkage_var),
    type="apeglm")
  resOrdered <- res[order(res$padj),]
} else {
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
}

# Write output to new file
write.table(as.data.frame(resOrdered), file = OUTPUT, sep = '\t', col.names = T, row.names = T)
