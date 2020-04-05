#!/usr/bin/env Rscript
license <- function() {
  "
  riboWaltz v1.1.0
  MIT License

  Copyright (c) [2017] [LAboratory of Translational Architectomics]

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the \"Software\"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
  "
}

description <- function() {
  "
  From BAM files to lists of data tables or GRangesList objects.

  This function reads one or multiple BAM files converting them into data
  tables or GRanges objects, arranged in a list or a GRangesList, respectively.
  In both cases the list elements contain, for each read: i) the name of the
  corresponding reference sequence (i.e. of the transcript on which it aligns);
  ii) its leftmost and rightmost position with respect to the 1st nucleotide of
  the reference sequence; iii) its length (intended as the width of the
  reference sequence region covered by the RNA fragment, see parameter
  indel_threshold and Details); iv) the leftmost and rightmost
  position of the annotated CDS of the reference sequence (if any) with respect
  to its 1st nucleotide. Please note: start and stop codon positions for
  transcripts without annotated CDS are set to 0.
  "
}

library(data.table)

bamtolist <- function(bamfolder, annotation, transcript_align = TRUE,
                      name_samples = NULL, indel_threshold = 5,
                      refseq_sep = NULL, granges = FALSE) {
  names <- list.files(path = bamfolder, pattern = ".bam$")
  if (length(name_samples) == 0) {
    name_samples <- unlist(strsplit(names, ".bam"))
    names(name_samples) <- unlist(strsplit(names, ".bam"))
  } else {
    if (length(name_samples) > length(names)) {
      cat("\n")
      stop("length of name_samples greater than number of files\n\n")
    }
    if (length(name_samples) < length(names)) {
      cat("\n")
      stop("length of name_samples smaller than number of files\n\n")
    }
  }

  i <- 0
  sample_reads_list <- list()
  for (n in names) {
    i <- i + 1
    cat(sprintf("Reading %s\n", n))
    filename <- paste(bamfolder, n, sep = "/")
    dt <- as.data.table(GenomicAlignments::readGAlignments(filename))
    nreads <- nrow(dt)
    cat(sprintf("Input reads: %s M\n", format(round((nreads / 1000000), 3), nsmall = 3)))
    dt <- dt[, diff_width := qwidth - width
             ][abs(diff_width) <= indel_threshold]
    if(nreads != nrow(dt)){
      perc_nreads <- round(((nreads - nrow(dt)) / nreads) * 100, 3)
      if(perc_nreads >= 0.001){
        cat(sprintf("%s M  (%s %%) reads removed: exceeding indel_threshold.\n",
                    format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                    format(perc_nreads, nsmall = 3) ))
      } else {
        cat(sprintf("%s (< 0.001 %%) reads removed: exceeding indel_threshold.\n",
                    format(round((nreads - nrow(dt)), 3), nsmall = 3)))
      }
    } else {
      cat("Good! Number of indel below indel_threshold for all reads. No reads removed.\n")
    }

    dt <- dt[, .(seqnames, start, end, width, strand)]
    setnames(dt, c("transcript", "end5", "end3", "length", "strand"))

    if(length(refseq_sep) != 0){
      dt <- dt[, transcript := tstrsplit(transcript, refseq_sep, fixed = TRUE, keep = 1)]
    }

    nreads <- nrow(dt)
    dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
    if(nreads != nrow(dt)){
      if(nrow(dt) == 0){
        stop(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n\n",
                     format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                     format(round(((nreads - nrow(dt)) / nreads) * 100, 3), nsmall = 3) ))
      } else{
        cat(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n",
                    format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                    format(round(((nreads - nrow(dt)) / nreads) * 100, 3), nsmall = 3) ))
      }
    } else {
      cat("Great! All reads' reference transcript IDs were found in the annotation table. No reads removed.\n")
    }

    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(dt)
      dt <- dt[strand == "+"]

      if(nreads != nrow(dt)){
        perc_nreads <- round(((nreads - nrow(dt)) / nreads) * 100, 3)
        if(perc_nreads >= 0.001){
          cat(sprintf("%s M  (%s %%) reads removed: mapping on negative strand.\n",
                      format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                      format(perc_nreads, nsmall = 3) ))
        } else {
          cat(sprintf("%s (< 0.001 %%) reads removed: mapping on negative strand.\n",
                      format(round((nreads - nrow(dt)), 3), nsmall = 3)))
        }
      } else {
        cat("Cool! All reads mapping on positive strand. No reads removed.\n")
      }
    }

    dt[annotation, on = 'transcript', c("cds_start", "cds_stop") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[cds_start == 1 & cds_stop == 0, cds_start := 0]
    dt[, strand := NULL]

    nreads <- nrow(dt)
    cat(sprintf("Output reads: %s M\n", format(round((nreads / 1000000), 3), nsmall = 3)))

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

    sample_reads_list[[name_samples[i]]] <- dt
    cat("\n")
  }

  if (granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }

  return(sample_reads_list)
}

#' From BAM files to BED files.
#'
#' This function reads one or multiple BAM files converting them into BED files
#' that contain, for each read: i) the name of the corresponding reference
#' sequence (i.e. of the transcript on which it aligns); ii) its leftmost and
#' rightmost position with respect to the 1st nucleotide of the reference
#' sequence; iii) its length (intended as the width of the reference sequence
#' region covered by the RNA fragment. For further information about this choice
#' please refer to section \code{Details} of function \code{\link{bamtolist}});
#' iv) the strand on which it aligns. Please note: this function relies on the
#' \emph{bamtobed} utility of the BEDTools suite and can be only run on UNIX,
#' LINUX and Apple OS X operating systems. Moreover, to generate R data
#' structures containing reads information, the \code{\link{bedtolist}} must be
#' run on the resulting BED files. For these reasons the authors suggest the use
#' of \code{\link{bamtolist}}.
#'
#' @param bamfolder Character string specifying the path to the folder storing
#'   BAM files. Please note: the function looks for BAM files recursively
#'   starting from the specified folder.
#' @param bedfolder Character string specifying the path to the directory where
#'   BED files shuold be stored. If the specified folder doesn't exist, it is
#'   automatically created. If NULL (the default), BED files are stored in a new
#'   subfolder of the working directory, called \emph{bed}.
#' @examples
#' ## path_bam <- "path/to/BAM/files"
#' ## path_bed <- "path/to/output/directory"
#' ## bamtobed(bamfolder = path_bam, bedfolder = path_bed)
#' @export
bamtobed <- function(bamfolder, bedfolder = NULL) {
  if (length(bedfolder) == 0) {
    bedfolder <- paste(bamfolder, "/bed", sep = "")
  }
  if (!dir.exists(bedfolder)) {
    dir.create(bedfolder)
  }
  bamtobed <- paste("bam_folder=\"", bamfolder, "\" && out_folder=\"", bedfolder, "\" && bam_files=$(find $bam_folder -type \"f\" -name \"*.bam\" | sort) && for name in $bam_files; do outname=$(echo $name | rev | cut -d'/' -f 1 | cut -d '.' -f 2- | rev); printf \"from \\t\\t $name \\t\\t to \\t\\t $out_folder/$outname.bed\\n\"; bedtools bamtobed -i $name | cut -f1,2,3,6 | awk -F'\\t' '{OFS = \"\\t\"; $5=$4; $2=$2+1; $4=$3-$2+1; print}' > $out_folder/$outname.bed; done",
                    sep = "")
  system(bamtobed)
}

#' From BED files to lists of data tables or GRangesList objects.
#'
#' This function reads one or multiple BED files, as generated by
#' \code{\link{bamtobed}}, converting them into data tables or GRanges objects,
#' arranged in a list or a GRangesList, respectively. In both cases two columns
#' are attached to the original data containing, for each read, the leftmost
#' and rightmost position of the annotated CDS of the reference sequence (if
#' any) with respect to its 1st nucleotide. Please note: start and stop codon
#' positions for transcripts without annotated CDS are set to 0.
#'
#' @param bedfolder Character string specifying the path to the folder storing
#'   BED files as generated by \code{\link{bamtobed}}.
#' @param annotation Data table as generated by \code{\link{create_annotation}}.
#'   Please make sure the name of reference transcripts in the annotation data
#'   table match those in the BED files (see also \code{refseq_sep}).
#' @param transcript_align Logical value whether BED files in \code{bedfolder}
#'   come from a transcriptome alignment (intended as an alignment against
#'   reference transcript sequences, see \code{Details}). If TRUE (the default),
#'   reads mapping on the negative strand should not be present and, if any,
#'   they are automatically removed.
#' @param name_samples Named character string vector specifying the desired name
#'   for the output list elements. A character string for each BED file in
#'   \code{bedfolder} is required. Plase be careful to name each element of the
#'   vector after the correct corresponding BED file in \code{bedfolder},
#'   leaving their path and extension out. No specific order is required.
#'   Default is NULL i.e. list elements are named after the name of the BED
#'   files, leaving their path and extension out.
#' @param refseq_sep Character specifying the separator between reference
#'   sequences' name and additional information to discard, stored in the same
#'   field (see \code{Details}). All characters before the first occurrence of
#'   the specified separator are kept. Default is NULL i.e. no string splitting
#'   is performed.
#' @param granges Logical value whether to return a GRangesList object. Default
#'   is FALSE i.e. a list of data tables is returned instead (the required input
#'   for \code{\link{length_filter}}, \code{\link{psite}},
#'   \code{\link{psite_info}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}).
#' @details \strong{riboWaltz} only works for read alignments based on
#'   transcript coordinates. This choice is due to the main purpose of RiboSeq
#'   assays to study translational events through the isolation and sequencing
#'   of ribosome protected fragments. Most reads from RiboSeq are supposed to
#'   map on mRNAs and not on introns and intergenic regions. Nevertheless, BAM
#'   based on transcript coordinates can be generated in two ways: i) aligning
#'   directly against transcript sequences; ii) aligning against standard
#'   chromosome sequences, requiring the outputs to be translated in transcript
#'   coordinates. The first option can be easily handled by many aligners (e.g.
#'   Bowtie), given a reference FASTA file where each sequence represents a
#'   transcript, from the beginning of the 5' UTR to the end of the 3' UTR. The
#'   second procedure is based on reference FASTA files where each sequence
#'   represents a chromosome, usually coupled with comprehensive gene annotation
#'   files (GTF or GFF). The STAR aligner, with its option --quantMode
#'   TranscriptomeSAM (see Chapter 6 of its
#'   \href{http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf}{manual}),
#'    is an example of tool providing such a feature.
#'
#'   \code{refseq_sep} is intended to lighten the identifiers of the reference
#'   sequences included in the final data table or to modify them to match those
#'   in the annotation table. Many details about the reference sequence such as
#'   their version (usually dot-separated), their length, name variants,
#'   associated gene/transcript/protein names (usually pipe-separated) might
#'   indeed be stored in the FASTA file used for the alignment and automatically
#'   transferred in the BAM.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' ## path_bed <- "path/to/BED/files"
#' ## bedtolist(bedfolder = path_bed, annotation = mm81cdna)
#' @import data.table
#' @export
bedtolist <- function(bedfolder, annotation, transcript_align = TRUE,
                      name_samples = NULL, refseq_sep = FALSE,
                      granges = FALSE) {
  names <- list.files(path = bedfolder, pattern = ".bed$")
  if (length(name_samples) == 0) {
    name_samples <- unlist(strsplit(names, ".bed"))
    names(name_samples) <- unlist(strsplit(names, ".bed"))
  } else {
    if (length(name_samples) > length(names)) {
      cat("\n")
      stop("length of name_samples greater than number of files\n\n")
    }
    if (length(name_samples) < length(names)) {
      cat("\n")
      stop("length of name_samples smaller than number of files\n\n")
    }
  }

  i <- 0
  sample_reads_list <- list()
  for (n in names) {
    i <- i + 1
    cat(sprintf("Reading %s\n", n))
    filename <- paste(bedfolder, n, sep = "/")
    dt <- fread(filename, sep = "\t", header = FALSE)
    setnames(dt, c("transcript", "end5", "end3", "length", "strand"))

    if(length(refseq_sep) != 0){
      dt <- dt[, transcript := tstrsplit(transcript, refseq_sep, fixed = TRUE, keep = 1)]
    }

    nreads <- nrow(dt)
    cat(sprintf("Input reads: %s M\n", format(round((nreads / 1000000), 3), nsmall = 3)))
    dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
    if(nreads != nrow(dt)){
      if(nrow(dt) == 0){
        stop(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n\n",
                     format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                     format(round(((nreads - nrow(dt)) / nreads) * 100, 3), nsmall = 3) ))
      } else{
        cat(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n",
                    format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                    format(round(((nreads - nrow(dt)) / nreads) * 100, 3), nsmall = 3) ))
      }
    } else {
      cat("Great! All reads' reference transcript IDs were found in annotation table. No reads removed.\n")
    }

    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(dt)
      dt <- dt[strand == "+"]

      if(nreads != nrow(dt)){
        perc_nreads <- round(((nreads - nrow(dt)) / nreads) * 100, 3)
        if(perc_nreads >= 0.001){
          cat(sprintf("%s M  (%s %%) reads removed: mapping on negative strand.\n",
                      format(round((nreads - nrow(dt)) / 1000000, 3), nsmall = 3),
                      format(perc_nreads, nsmall = 3) ))
        } else {
          cat(sprintf("%s (< 0.001 %%) reads removed: mapping on negative strand.\n",
                      format(round((nreads - nrow(dt)), 3), nsmall = 3)))
        }
      } else {
        cat("Cool! All reads mapping on positive strand. No reads removed.\n")
      }
    }

    dt[annotation, on = 'transcript', c("cds_start", "cds_stop") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[cds_start == 1 & cds_stop == 0, cds_start := 0]
    dt[, strand := NULL]

    nreads <- nrow(dt)
    cat(sprintf("Output reads: %s M\n", format(round((nreads / 1000000), 3), nsmall = 3)))

    if(granges == T || granges == TRUE) {
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

    sample_reads_list[[name_samples[i]]] <- dt
  }

  if(granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }

  return(sample_reads_list)
}
