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

# Import dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org")}

if ("data.table" %in% rownames(installed.packages()) == FALSE) {
  print("Installing data.table...")
  install.packages("data.table", repos = "http://cran.us.r-project.org")
} else {
  print("data.table package already installed")
}
library(data.table)

if ("GenomicAlignments" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("GenomicAlignments", dependencies=c("Depends", "Imports", "LinkingTo"))
} else {
  print("GenomicAlignments package already installed")
}
library(GenomicAlignments)

# Set globals
chromosome <- 'chromosome'
left_coordinate <- 'left_coordinate'
right_coordinate <- 'right_coordinate'
strand <- 'strand'
exon <- 'exon'
sample_factor <- 0.1

# Get arguments
# args[1] = Path to RPF BAM file
# args[2] = List of Cufflinks isoform counts table
# args[3] = Output file path
args = commandArgs(trailingOnly=TRUE)

# Set parameters
RPF_FILE <- args[1]
RNA_FILE <- args[2]
OUTPUT_LOCATION <- args[3]

# Parse cufflinks file

# Make transcript dictionary {gene_id: most_abundant_transcript_id}

# Init count table for isoforms based on cufflinks file

# Parse out any seq records not in frame?

# Count up unique aligners

# Collapse multimappers based on most abundant and add counts to table

# Output table to file
