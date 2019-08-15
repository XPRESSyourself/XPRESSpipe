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

# Measure library complexity of RNA-seq sample

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

if ("Rsubread" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("Rsubread", dependencies=c("Depends", "Imports", "LinkingTo"))
} else {
  print("Rsubread package already installed")
}
if ("dupRadar" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("dupRadar", dependencies=c("Depends", "Imports", "LinkingTo"))
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
