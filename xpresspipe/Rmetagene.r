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
sample_factor <- 0.1

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

      # Get representative sample of BAM file
      rows <- dim(bam)[1]
      sample_amount <- ceiling(rows * sample_factor)
      bam <- bam[sample(nrow(bam), sample_amount), ]

      print(paste('Processing', toString(sample_amount), 'BAM records', toString(sample_factor), '(percent of total reads)...', sep=' '))
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
process_metacoverage <- function(
  bam, index) {

    # Get records lengths and find unique transcripts to parse through
    index$length <- index$right_coordinate - index$left_coordinate + 1
    index_exon <- index[index$feature == exon,]
    unique_transcripts <- as.character(unique(unlist(index_exon$transcript)))

    # Get total length of exon space
    transcript_dictionary <- data.frame(
      'number' = 1:length(unique_transcripts),
      'transcript' = unique_transcripts,
      'total_exon_length' = '',
      row.names = 'number',
      stringsAsFactors = FALSE)

    for (transcript_id in unique_transcripts) {

      # Calculate exon space for the gene
      cumul_exon_length <- sum(index_exon[index_exon$transcript == transcript_id,]$length)

      if (total_exon_length == 0) {

      } else {

          # Add the length to the dictionary
          transcript_dictionary[which(transcript_dictionary$name == transcript_id),'total_exon_length'] <- cumul_exon_length
      }
    }

    # Make empty dataframe with min/max range
    metacounts <- data.frame('position' = 1:100, 'metacount' = 0, row.names = 'position', stringsAsFactors = FALSE)

    # Get meta-positions from BAM file
    bam_filter <- bam[bam$seqnames %in% unique_transcripts]
    bam_filter$meta_position <- bam_filter$end - floor((bam_filter$end - bam_filter$start) / 2)

    # Loop through start and end of each read in range and add one for every nt position
    for (read in 1:nrow(bam_filter)) {

      if (bam_filter[read,'seqnames'] %in% unique_transcripts) {

        # Add a count for each metaposition
        exon_len <- transcript_dictionary[which(transcript_dictionary$name == bam_filter[read,'seqnames']),'total_exon_length']
        meta_pos <- ceiling(bam_filter[read,'meta_position'] / exon_len)
        metacounts[meta_pos, 'metacount'] <- metacounts[meta_pos, 'metacount'] + 1

      } else {
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

      # Analyze
      if (dim(bam)[1] == 0) {

        print(paste('No records for', toString(target_transcript), 'were found', sep=' '))

      } else {

        coverage <- process_metacoverage(bam, index)
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
