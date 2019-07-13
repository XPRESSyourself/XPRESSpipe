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

if ("doMC" %in% rownames(installed.packages()) == FALSE) {
  print("Installing doMC...")
  install.packages("doMC", repos = "http://cran.us.r-project.org")
} else {
  print("doMC package already installed")
}
library(doMC)

library(parallel)


# Set globals
chromosome <- 'chromosome'
left_coordinate <- 'left_coordinate'
right_coordinate <- 'right_coordinate'
strand <- 'strand'
buffer <- 1000
n_cores <- detectCores()

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

      print('Splitting BAM record into chunks...')
      rows <- nrow(dt)
      chunks <- split(
        bam,
        rep(1:ceiling(rows/n_cores),
        each=n_cores,
        length.out=rows))

      print(paste('Processing', rows, 'BAM records...', sep=' '))
      rm(dt)
      gc()
      return(chunks)

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

    # Get BAM file within range
    min <- min(index$left_coordinate) - buffer
    max <- max(index$right_coordinate) + buffer
    gene_strand <- levels(droplevels(index$strand))[1]
    chr <- index$chromosome[1]

    # Get coverage regions
    coverage_ranges <- character(0)
    for (record in 1:nrow(index)) {
      coverage_ranges <- c(coverage_ranges, toString(index[record,left_coordinate]):toString(index[record,right_coordinate]))
    }

    # Get relevant data
    bam_sub <- bam[(bam$seqnames == chr) & (bam$start >= min) & (bam$end <= max) & (bam$strand == gene_strand)]

    # Get overlapping coverage in region
    # Make empty dataframe with min/max range
    counts <- data.frame('position' = toString(min):toString(max), 'coverage' = 0, row.names = 'position')

    # Loop through start and end of each read in range and add one for every nt position
    for (read in 1:nrow(bam_sub)) {

      # Add a count for each nucleotide of the read
      counts[row.names(counts) %in% toString(bam_sub[read,start]):toString(bam_sub[read,end]),] <- (
        counts[row.names(counts) %in% toString(bam_sub[read,start]):toString(bam_sub[read,end]),] + 1
        )

      }

    # Keep only nt indices that fall in an exon range from index
    # Make list of intron exon positions
    # Stay on the same gene name
    # Only keep those indices from counts
    coverage_counts <- data.frame(
      'position' = coverage_ranges,
      'coverage' = counts[row.names(counts) %in% coverage_ranges,],
      row.names = 'position')

    return(coverage_counts)

  }

run_list <- function (
  file_path, file_list, index_file, output_path) {

    index <- fetch_index(index_file)

    for (f in file_list) {

      # Import BAM and get coverage
      file <- paste(file_path, f, sep='')
      chunks <- read_bam(file)

      # Multi-process BAM chunks
      registerDoMC()
      coverage <- foreach(r = chunks) %dopar% process_coverage(r, index)

      # Write BAM coverage metrics to output file
      file_name = vapply(strsplit(f, "[.]"), `[`, 1, FUN.VALUE=character(1))
      output_file = paste(output_path, file_name, '_metrics.txt', sep='')
      write.table(coverage, file=output_file, sep='\t', col.names=NA)

      # Clean the batch
      rm(file)
      rm(bam)
      rm(coverage)
      rm(file_name)
      rm(output_file)
      gc()

    }

  }

# MAIN
# Run files through coverage checker
run_list(BAM_PATH, BAM_LIST, INDEX, OUTPUT_LOCATION)
