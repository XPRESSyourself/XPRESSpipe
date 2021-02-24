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
install.packages("stringi", repos = "http://cran.us.r-project.org", quiet = TRUE)
library(stringi)

if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)}

if ("data.table" %in% rownames(installed.packages()) == FALSE) {
  print("Installing data.table...")
  install.packages("data.table", repos = "http://cran.us.r-project.org", quiet = TRUE)
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
sample_factor <- 1

# Get arguments
# args[1] = Path to BAM files
# args[2] = List of BAM files
# args[4] = Output file path
args = commandArgs(trailingOnly=TRUE)

# Set parameters
BAM_PATH <- args[1]
BAM_LIST <- strsplit(args[2],',')[[1]]
OUTPUT_LOCATION <- args[3]

# func: Import BAM file
# @param bam_file: BAM-format file
# @return: GenomicAlignments data.frame
read_bam <- function(
  bam_file) {

    if (endsWith(bam_file, '.bam')) {

      # Read in file to Genomic Alignments data.frame and report number of records
      print(paste('Reading', bam_file, sep=' '))
      bam <- as.data.table(GenomicAlignments::readGAlignments(bam_file))

      # Get representative sample of BAM file
      rows <- dim(bam)[1]
      sample_amount <- ceiling(rows * sample_factor)
      bam <- bam[sample(nrow(bam), sample_amount), ]

      print(paste('Processing', toString(sample_amount), 'BAM records (', toString(sample_factor * 100), 'percent of total reads )...', sep=' '))
      return(bam)

    } else {

      # Raise warning if file is not BAM file
      print(c(bam_file, ' is not a BAM file, skipping...'))
    }
  }

# Requires a STAR styled transcriptome-aligned BAM file
# Note: STAR transcriptome-aligned BAM files only give alignments that are within processed mRNA space, so introns are not included by default, therefore coordinates fall within UTRs or CDSs
# https://www.biorxiv.org/content/biorxiv/early/2018/10/16/444620.full.pdf
run_list <- function (
  file_path, file_list, output_path) {

    for (f in file_list) {

      # Import BAM and get coverage
      print(paste('Processing records for', toString(f), sep=' '))
      file <- paste(file_path, f, sep='')
      print(paste('Reading records for', toString(f), sep=' '))
      bam <- read_bam(file)

      # Analyze
      if (dim(bam)[1] == 0) {

        print(paste('No records for', toString(target_transcript), 'were found', sep=' '))

      } else {

        bam$meta_position <- bam$end - floor((bam$end - bam$start) / 2)
        bam <- bam[,c('seqnames','meta_position')]

      }

      # Write BAM coverage metrics to output file
      file_name = vapply(strsplit(f, "[.]"), `[`, 1, FUN.VALUE=character(1))
      output_file = paste(output_path, file_name, '.metaposit', sep='')
      write.table(bam, file=output_file, sep='\t', col.names=NA)

      print(paste('Metagene profile for', toString(f), 'written to', toString(output_file), sep=' '))
    }
  }

# MAIN
# Run files through coverage checker
print(paste('Rscript processing file list:', BAM_LIST, sep=' '))
run_list(BAM_PATH, BAM_LIST, OUTPUT_LOCATION)
print("Metagene BAM file processing complete.")
