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
if ("data.table" %in% rownames(installed.packages()) == FALSE) {
  print("Installing data.table...")
  install.packages("data.table", repos = "http://cran.us.r-project.org")
} else {
  print("data.table package already installed")
}
library(data.table)

if ("GenomicAlignments" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install("GenomicAlignments", dependencies=TRUE)
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

# func: Get appropriate range of BAM file based on index and count coverage
# @param bam: GenomicAlignments data.frame
# @param index: XPRESSpipe-formatted index data.frame object
# @return: Coverage data.frame
process_coverage <- function(
  bam, index) {

    # Normalize coordinate positions to match BAM file, start position starts at 1
    min_coordinate <- min(index$left_coordinate) - 1
    index$right_coordinate <- index$right_coordinate - min_coordinate
    index$left_coordinate <- index$left_coordinate - min_coordinate
    index$length <- index$right_coordinate - index$left_coordinate + 1

    # Get total length of exon space
    index_exon <- index[index$feature == exon,]
    transcript_length <- sum(index_exon$length)

    # Get overlapping coverage in region
    # Make empty dataframe with min/max range
    counts <- data.frame('position' = 1:transcript_length, 'coverage' = 0, 'feature' = '', row.names = 'position', stringsAsFactors = FALSE)

    # Loop through start and end of each read in range and add one for every nt position
    for (read in 1:nrow(bam)) {

      # Add a count for each nucleotide of the read if part of the transcript (will allow partial mapping)
      counts[row.names(counts) %in% toString(bam[read,start]):toString(bam[read,end]),'coverage'] <- (
        counts[row.names(counts) %in% toString(bam[read,start]):toString(bam[read,end]),'coverage'] + 1
        )

      }

    # Get exon start coordinates
    for (record in 1:nrow(index_exon)) {
      if (record == 1) {
        index_exon[record,'exon_start'] <- 1
      } else {

        index_exon[record,'exon_start'] <- index_exon[record - 1,'length'] + index_exon[record - 1,'exon_start']

      }

    }

    # Label exon feature starts
    exon_number <- 1
    for (record in 1:nrow(counts)) {

      if (record %in% index_exon$exon_start) {

        # Add exon label to start of each exon iteratively
        counts[record,'feature'] <- paste('Exon', toString(exon_number), sep=' ')
        exon_number <- exon_number + 1

      }

    }

    # Keep only nt indices that fall in an exon range from index
    # Make list of intron exon positions
    # Stay on the same gene name
    # Only keep those indices from counts
    #coverage_counts <- data.frame(
    #  'position' = 1:nrow(counts),
    #  'coverage' = counts$coverage,
    #  'feature' = counts$feature,
    #  row.names = 'position')

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

        coverage <- process_coverage(target_bam, index)

      }
      rm(file)
      rm(bam)
      gc()

      # Write BAM coverage metrics to output file
      file_name = vapply(strsplit(f, "[.]"), `[`, 1, FUN.VALUE=character(1))
      output_file = paste(output_path, file_name, '_metrics.txt', sep='')
      write.table(coverage, file=output_file, sep='\t', col.names=NA)

      # Clean the batch
      rm(target_bam)
      rm(coverage)
      rm(file_name)
      rm(output_file)
      gc()

    }

  }

# MAIN
# Run files through coverage checker
run_list(BAM_PATH, BAM_LIST, INDEX, OUTPUT_LOCATION)
