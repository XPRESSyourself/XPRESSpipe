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
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("BH", quietly = TRUE)) {install.packages("BH", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("cpp11", quietly = TRUE)) {install.packages("cpp11", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("plogr", quietly = TRUE)) {install.packages("plogr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("XML", quietly = TRUE)) {install.packages("XML", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("dbplyr", quietly = TRUE)) {install.packages("dbplyr", repos = "http://cran.us.r-project.org", quiet = TRUE)}
if (!requireNamespace("data.table", quietly = TRUE)) {install.packages("data.table", repos = "http://cran.us.r-project.org", quiet = TRUE)}


if ("GenomicAlignments" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("GenomicAlignments", dependencies=c("Depends", "Imports"), quiet = TRUE)
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

# Get arguments
# args[1] = Path to BAM files
# args[2] = List of BAM files
# args[3] = Index file with path
# args[4] = Output file path
args = commandArgs(trailingOnly=TRUE)

# Set parameters
BAM_PATH <- args[1]
BAM_LIST <- strsplit(args[2],',')[[1]]
INDEX <- args[3]
OUTPUT_LOCATION <- args[4]

# func: Import BAM file
# @param bam_file: BAM-format file
# @return: GenomicAlignments data.frame
read_bam <- function(
  bam_file) {

    if (endsWith(bam_file, '.bam')) {

      # Read in file to Genomic Alignments data.frame and report number of records
      print(paste('Reading', bam_file, sep=' '))
      bam <- as.data.table(GenomicAlignments::readGAlignments(bam_file))
      rows <- dim(bam)[1]
      print(paste('Processing', rows, 'BAM records...', sep=' '))
      return(bam)

    } else {

      # Raise warning if file is not BAM file
      print(c(bam_file, ' is not a BAM file, skipping...'))

    }

  }

# func: Import index file
# @param index_file: XPRESSpipe-formatted index file
# @return: Index data.frame
fetch_index <- function(
  index_file) {

    # Read in gene index
    index <- read.table(
      index_file,
      header = TRUE,
      sep = '\t')

    return(index)

  }

# func: Get appropriate range of BAM file based on index and count genecoverage
# @param bam: GenomicAlignments data.frame
# @param index: XPRESSpipe-formatted index data.frame object
# @return: geneCoverage data.frame
# Note: STAR transcriptome-aligned BAM files only give alignments that are within processed mRNA space, so introns are not included by default, therefore coordinates fall within UTRs or CDSs
# https://www.biorxiv.org/content/biorxiv/early/2018/10/16/444620.full.pdf
process_coverage <- function(
  bam, index) {

    # Get overlapping coverage in region
    # Make empty dataframe with min/max range
    transcript_length <- sum(index$l_tr)
    counts <- data.frame('position' = 1:transcript_length, 'coverage' = 0, row.names = 'position', stringsAsFactors = FALSE)

    # Loop through start and end of each read in range and add one for every nt position
    for (read in 1:nrow(bam)) {

      # Add a count for each nucleotide of the read if part of the transcript (will allow partial mapping)
      counts[row.names(counts) %in% toString(bam[read,start]):toString(bam[read,end]),'coverage'] <- (
        counts[row.names(counts) %in% toString(bam[read,start]):toString(bam[read,end]),'coverage'] + 1
        )

      }

    return(counts)

  }

# Requires a STAR styled transcriptome-aligned BAM file
run_list <- function (
  file_path, file_list, index_file, output_path) {

    index <- fetch_index(index_file)

    for (f in file_list) {

      # Import BAM and get coverage
      file <- paste(file_path, f, sep='')
      bam <- read_bam(file)

      # Get transcript to analyze
      target_transcript <- unique(levels(droplevels(index$transcript)))
      target_bam <- bam[bam$seqnames == target_transcript]
      if (dim(target_bam)[1] == 0) {

        print(paste('No records for', toString(target_transcript), 'were found', sep=' '))

      } else {

        genecoverage <- process_coverage(target_bam, index)

      }
      rm(file)
      rm(bam)
      gc()

      # Write BAM coverage metrics to output file
      file_name = vapply(strsplit(f, "[.]"), `[`, 1, FUN.VALUE=character(1))
      output_file = paste(output_path, file_name, '_metrics.txt', sep='')
      write.table(genecoverage, file=output_file, sep='\t', col.names=NA)

      # Clean the batch
      rm(target_bam)
      rm(genecoverage)
      rm(file_name)
      rm(output_file)
      gc()

    }

  }

# MAIN
# Run files through coverage checker
run_list(BAM_PATH, BAM_LIST, INDEX, OUTPUT_LOCATION)
