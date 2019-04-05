#!/usr/bin/env Rscript

# Measure periodicity of footprint sample

# Install dependencies
install.packages("devtools")
library(devtools)
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE,
				build_opts = c("--no-resave-data", "--no-manual"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicAlignments", version = "3.8", dependencies = TRUE)
BiocManager::install("GenomicRanges", version = "3.8", dependencies = TRUE)
BiocManager::install("GenomicFeatures", version = "3.8", dependencies = TRUE)
BiocManager::install("Rsamtools", version = "3.8", dependencies = TRUE)
BiocManager::install("ensembldb", version = "3.8")
library(riboWaltz)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(ensembldb)

# Get arguments
#
#
# args[1] = if TRUE, make and output reference
# args[2] = reference input --
# 	if args[1] is True, path and filename of transcripts.gtf
# 	if args[1] is False, path and filename of transcripts_flat.txt
# args[3] = path to bam files to process
# args[4] = output directory for metrics
args = commandArgs(trailingOnly=TRUE)

# Make reference for periodicity plotting
make_reference <- function(reference_file) {

	# Create output file name
	output <- paste(substr(reference_file, 0, (nchar(args[6])-4)), '_flat.txt', sep='');

	# Generate reference annotation for periodicity plotting
	gtf_file <- create_reference(reference_file);

	# Output annotation table
	write.table(as.data.frame(gtf_file), file=output, sep='\t', col.names=T, row.names=T);
}

# Prepare bam file and annotation for analysis
# Adapted from original riboWaltz code to allow for singular file input
bam_prep <- function(bamfile, annotation, transcript_align = TRUE,
											indel_threshold = 5,
											refseq_sep = NULL, granges = FALSE) {

	# Process single BAM file
	sample_reads_list <- list() # Initialize output list for feeding into normal riboWaltz functions
	dt <- as.data.frame(GenomicAlignments::readGAlignments(bamfile)) # Read BAM file into dataframe
	nreads <- nrow(dt) # Get number of records in

	dt[, "diff_width"] <- (dt[, "qwidth"] - dt[, "width"]) # Find difference of values
	dt[, "diff_width"] <- abs(dt[, "diff_width"]) # Convert diff_width to absolute values
	dt <- dt[dt$diff_width <= indel_threshold, ] # Threshold diff_width based off indel_threshold

	if(nreads != nrow(dt)){ # Remove reads not meeting threshold criteria
		cat(sprintf("%s M  (%s %%) reads removed: exceeding indel_threshold.\n",
               	format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2),
               	format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
	} else {
		cat("Number of indel below indel_threshold for all reads. No reads removed.\n")
 	}

 	dt <- dt[, c("seqnames", "start", "end", "width", "strand")]
 	names(dt) <- c("transcript", "end5", "end3", "length", "strand")



 	nreads <- nrow(dt)
 	cat(sprintf("Reads: %s M\n", format(round((nreads / 1000000), 2), nsmall = 2)))
 	dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
 	if(nreads != nrow(dt)){
   	if(nrow(dt) == 0){
     	stop(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n\n",
                  	format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2),
                  	format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
   	} else{
    	cat(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n",
                 	format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2),
                 	format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
   	}
 	} else {
   	cat("All reads' reference transcript IDs were found in annotation table. No reads removed.\n")
 	}

 	if(transcript_align == TRUE | transcript_align == T){
   	nreads <- nrow(dt)
   	dt <- dt[strand == "+"]

   	if(nreads != nrow(dt)){
     	cat(sprintf("%s M  (%s %%) reads removed: mapping on negative strand.\n",
                 	format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2),
                 	format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
   	} else {
     	cat("All reads mapping on positive strand. No reads removed.\n")
   	}
 	}

 	dt[annotation, on = 'transcript', c("cds_start", "cds_stop") <- list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
 	dt[cds_start == 1 & cds_stop == 0, cds_start <- 0]
 	dt[, strand <- NULL]

 	if (granges == T || granges == TRUE) {
   	dt <- GenomicRanges::makeGRangesFromDataFrame(dt,
                                                 	keep.extra.columns = TRUE,
                                                 	ignore.strand = TRUE,
                                                 	seqnames.field = c("transcript"),
                                                 	start.field = "end5",
                                                 	end.field = "end3",
                                                 	strand.field = "strand",
                                                 	starts.in.df.are.0based = FALSE)
   	GenomicRanges::strand(dt) <- "+"
 	}

 	sample_reads_list[[sampname]] <- dt

 	if (granges == T || granges == TRUE) {
   	sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
 	}

 	return(sample_reads_list)
}

# Generate periodicity metrics
make_periodicity <- function(reference_annotation, bam_file, output) {

	# Read in annotation file
	annotation_dt <- read.table(file = reference_annotation, sep = '\t', header = TRUE)

	# Read in BAM file to process
	bam_list <- bam_prep(bamfile = bam_file, annotation = annotation_dt)

	filtered_list <- length_filter(data = bam_list, length_filter_mode = "custom",
				 length_filter_vector = 27:30)

}

# MAIN::Make reference for riboWaltz if arguments satisfied
if (isTRUE(args[1]) {
	make_reference(args[2])
}
else {
	make_periodicity()
}
