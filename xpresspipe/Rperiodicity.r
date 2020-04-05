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

# Run install of data.table to make sure most recent version is installed
library(data.table)

# Get arguments
# args[1] = BAMtranscriptome directory
# args[2] = GTF reference file (full)
# args[3] = output directory for metrics
# args[4] = file directory
args = commandArgs(trailingOnly=TRUE)

# Set parameters
BAM_LIST <- args[1] # Directory containing
GTF <- args[2]
OUTPUT_P_SITES <- args[3] # Path and filename with .txt extension
PATH <- args[4] # File directory

# Install riboWaltz functions
source(paste(toString(PATH), "ribowaltz_annotation.R", sep=""))
source(paste(toString(PATH), "ribowaltz_psites.R", sep=""))
source(paste(toString(PATH), "ribowaltz_combine_bam.R", sep=""))

# Get p-site offsets
annotation_dt <- create_annotation(gtfpath=GTF)

# Get data for codon usage -- output to both directories
write.table(
  as.data.frame(annotation_dt),
  file=paste(
    OUTPUT_P_SITES,
    "annotation.txt",
    sep=""),
  sep='\t',
  col.names=NA)

# Get p-sites
reads_list <- bamtolist(bamfolder = BAM_LIST, annotation = annotation_dt)
p_sites <- psite(reads_list, extremity="5end") # This will fail if the input files are too small and don't have good representation across genes

if (typeof(p_sites) != "NULL") {
  p_info <- psite_info(reads_list, p_sites)

  # Get list of unique elements in 'sample' column in p_sites
  # Generate tables for each sample
  for (SAMPLE in as.list(unique(p_sites$sample))) {
    SAMPLE_NAME = vapply(strsplit(SAMPLE, "[.]"), `[`, 1, FUN.VALUE=character(1)) # Get sample name

    OUTPUT_NAME = paste(OUTPUT_P_SITES, SAMPLE_NAME, "_metrics.txt", sep="")

    write.table(as.data.frame(
      p_info[[SAMPLE]][,c("transcript","psite","length")]
    ), file=OUTPUT_NAME, sep='\t', col.names=NA)
  }
}
