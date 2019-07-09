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
# args[3] = output directory for metrics
args = commandArgs(trailingOnly=TRUE)

# Set parameters
BAM_LIST <- args[1] # Directory containing
GTF <- args[2]
OUTPUT <- args[3] # Path and filename with .txt extension
print(args[1])
print(args[2])
print(args[3])


# Get p-site offsets
annotation_dt <- create_annotation(gtfpath = GTF)
reads_list <- bamtolist(bamfolder = BAM_LIST, annotation = annotation_dt)
p_sites <- psite(reads_list)
p_info <- psite_info(reads_list, p_sites)

# Get list of unique elements in 'sample' column in p_sites
# Generate tables for each sample
for (SAMPLE in as.list(unique(p_sites$sample))) {
  SAMPLE_NAME = vapply(strsplit(SAMPLE, "[.]"), `[`, 1, FUN.VALUE=character(1))
  OUTPUT_NAME = paste(OUTPUT, SAMPLE_NAME, "_metrics.txt", sep="")
  meta_prof <- metaprofile_psite(p_info, annotation_dt, SAMPLE, cdsl = 75)
  write.table(as.data.frame(meta_prof$dt), file=OUTPUT_NAME, sep='\t', col.names=NA)
}
