#!/usr/bin/env Rscript

# Measure library complexity of RNA-seq sample

# Install dependencies
if ("devtools" %in% rownames(installed.packages()) == FALSE) {
  print("Installing devtools...")
  install.packages("devtools", repos = "http://cran.us.r-project.org")
} else {
  print("devtools package already installed")
}
library(devtools)

if ("riboWaltz" %in% rownames(installed.packages()) == FALSE) {
  print("Installing riboWaltz...")
  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
} else {
  print("riboWaltz package already installed")
}
library(riboWaltz)

# Get arguments
# args[1] = BAMtranscriptome directory
# args[2] = GTF reference file (full)
# args[3] = output file for metrics
args = commandArgs(trailingOnly=TRUE)

# Set parameters
bams <- args[1]
gtf <- args[2]
output <- args[3] # Path and filename with .txt extension

annotation_dt <- create_annotation(gtfpath = gtf)
reads_list <- bamtolist(bamfolder = bams, annotation = annotation_dt)
p_sites <- psite(reads_list)

write.table(as.data.frame(p_sites), file=output, sep='\t', col.names=NA)
