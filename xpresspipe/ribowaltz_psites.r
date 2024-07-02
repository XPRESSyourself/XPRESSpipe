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
  Ribosome P-sites position within reads.

  This function identifies the exact position of the ribosome P-site within
  each read, determined by the localisation of its first nucleotide (see
  Details). It returns a data table containing, for all samples and read
  lengths: i) the percentage of reads in the whole dataset, ii) the percentage
  of reads aligning on the start codon (if any); iii) the distance of the
  P-site from the two extremities of the reads before and after the correction
  step; iv) the name of the sample. Optionally, this function plots a
  collection of read length-specific occupancy metaprofiles displaying the
  P-site offsets computed through the process.
  "
}

library(data.table)

psite <- function(
    data,
    flanking=6,
    start=TRUE,
    extremity="auto",
    plot=FALSE,
    plot_dir=NULL,
    plot_format="png",
    cl=99,
    log_file=FALSE,
    log_file_dir=NULL) {

  if (log_file) {
    if (length(log_file_dir) == 0) {
      log_file_dir <- getwd()
    }
    if (!dir.exists(log_file_dir)) {
      dir.create(log_file_dir)
    }
    logpath <- paste0(log_file_dir, "/best_offset.txt")
    cat("sample\texremity\toffset(nts)\n", file = logpath)
  }

  names <- names(data)
  offset <- NULL
  
  for (n in names) {
    dt <- data[[n]]
    
    # Ensure necessary columns are present
    required_columns <- c("transcript", "end5", "end3", "cds_start", "cds_stop", "length", "l_utr5", "l_utr3")
    missing_columns <- setdiff(required_columns, colnames(dt))
    
    if (length(missing_columns) > 0) {
      stop(sprintf("Missing required columns: %s", paste(missing_columns, collapse = ", ")))
    }

    lev <- sort(unique(dt$length))
    
    if (start) {
      base <- 0
      dt[, site_dist_end5 := ifelse(is.na(cds_start), NA, end5 - cds_start)]
      dt[, site_dist_end3 := ifelse(is.na(cds_start), NA, end3 - cds_start)]
    } else {
      base <- -5
      dt[, site_dist_end5 := ifelse(is.na(cds_stop), NA, end5 - cds_stop - base)]
      dt[, site_dist_end3 := ifelse(is.na(cds_stop), NA, end3 - cds_stop - base)]
    }
    
    site_sub <- dt  #[!is.na(site_dist_end5) & site_dist_end5 <= -flanking & !is.na(site_dist_end3) & site_dist_end3 >= flanking - 1]


    if (nrow(site_sub) == 0) {
      warning(sprintf("Could not find valid records from %s", n))
      next
    }
    
    minlen <- min(site_sub$length, na.rm = TRUE)
    maxlen <- max(site_sub$length, na.rm = TRUE)

    t <- table(factor(site_sub$length, levels = lev))

    offset_temp <- data.table(length = as.numeric(as.character(names(t))), percentage = (as.vector(t) / sum(as.vector(t))) * 100)
    offset_temp[, around_site := "T"
                ][percentage == 0, around_site := "F"]
    
    tempoff <- function(v_dist) {
      ttable <- sort(table(v_dist), decreasing = TRUE)
      ttable_sr <- ttable[as.character(as.numeric(names(ttable)) + 1)]
      ttable_sl <- ttable[as.character(as.numeric(names(ttable)) - 1)]
      tsel <- rowSums(cbind(ttable > ttable_sr, ttable > ttable_sl), na.rm = TRUE)
      return(as.numeric(names(tsel[tsel == 2][1])))
    }

    offset_temp5 <- site_sub[, list(offset_from_5 = tempoff(.SD$site_dist_end5)), by = length]
    offset_temp3 <- site_sub[, list(offset_from_3 = tempoff(.SD$site_dist_end3)), by = length]
    merge_allx <- function(x, y) merge(x, y, all.x = TRUE, by = "length")
    offset_temp <- Reduce(merge_allx, list(offset_temp, offset_temp5, offset_temp3))

    adj_off <- function(dt_site, dist_site, add, bestoff) {
      temp_v <- dt_site[[dist_site]]
      t <- table(factor(temp_v, levels = seq(min(temp_v, na.rm = TRUE) - 2, max(temp_v, na.rm = TRUE) + add)))
      t[1:2] <- t[3] + 1
      locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
      adjoff <- locmax[which.min(abs(locmax - bestoff))]
      ifelse(length(adjoff) != 0, adjoff, bestoff)
    }

    best_from5_tab <- offset_temp[!is.na(offset_from_5), list(perc = sum(percentage)), by = offset_from_5
                                  ][perc == max(perc)]
    best_from3_tab <- offset_temp[!is.na(offset_from_3), list(perc = sum(percentage)), by = offset_from_3
                                  ][perc == max(perc)]

    if ((nrow(best_from3_tab) == 0) | (nrow(best_from5_tab) == 0)) {
      warning(sprintf("Could not find valid records from %s", n))
    } else {
      if (extremity == "auto" &
         ((best_from3_tab[, perc] > best_from5_tab[, perc] &
           as.numeric(best_from3_tab[, offset_from_3]) <= minlen - 2) |
          (best_from3_tab[, perc] <= best_from5_tab[, perc] &
           as.numeric(best_from5_tab[, offset_from_5]) > minlen - 1)) |
         extremity == "3end") {
        best_offset <- as.numeric(best_from3_tab[, offset_from_3])
        line_plot <- "3end"
        adj_tab <- site_sub[, list(corrected_offset_from_3 = adj_off(.SD, "site_dist_end3", 0, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
        offset_temp[is.na(corrected_offset_from_3), corrected_offset_from_3 := best_offset
                    ][, corrected_offset_from_5 := -corrected_offset_from_3 + length - 1]
      } else {
        if(extremity == "auto" &
           ((best_from3_tab[, perc] <= best_from5_tab[, perc] &
             as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1) |
            (best_from3_tab[, perc] > best_from5_tab[, perc] &
             as.numeric(best_from3_tab[, offset_from_3]) > minlen - 2)) |
           extremity == "5end"){
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "5end"
        adj_tab <- site_sub[, list(corrected_offset_from_5 = adj_off(.SD, "site_dist_end5", 1, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
        offset_temp[is.na(corrected_offset_from_5), corrected_offset_from_5 := best_offset
                    ][, corrected_offset_from_5 := abs(corrected_offset_from_5)
                      ][, corrected_offset_from_3 := abs(corrected_offset_from_5 - length + 1)]
        }
      }

      cat(sprintf("best offset: %i nts from the %s\n", abs(best_offset), gsub("end", "' end", line_plot)))

      if (log_file) {
        cat(sprintf("%s\t%s\t%i\n", n, gsub("end", "'end", line_plot), abs(best_offset)), file = logpath, append = TRUE)
      }

      t <- table(factor(dt$length, levels = lev))
      offset_temp[!is.na(offset_from_5), offset_from_5 := abs(offset_from_5)
                  ][, total_percentage := as.numeric(format(round((as.vector(t) / sum(as.vector(t))) * 100, 3), nsmall = 4))
                    ][, percentage := as.numeric(format(round(percentage, 3), nsmall = 4))
                      ][, sample := n]

      setcolorder(offset_temp, c("length", "total_percentage", "percentage", "around_site", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
      if (start) {
        setnames(offset_temp, c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
      } else {
        setnames(offset_temp, c("length", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
      }

      dt[, c("site_dist_end5", "site_dist_end3") := NULL]
      offset <- rbind(offset, offset_temp)
    }
  }
  return(offset)
}

description2 <- function() {
  "
  #' Update reads information according to the inferred P-sites.
  #'
  #' This function provides additional reads information according to the position
  #' of the P-site identfied by psite. It attaches to each data
  #' table in a list four columns reporting i) the P-site position with respect to
  #' the 1st nucleotide of the transcript, ii) the P-site position with respect to
  #' the start and the stop codon of the annotated coding sequence (if any) and
  #' iii) the region of the transcript (5' UTR, CDS, 3' UTR) that includes the
  #' P-site. Please note: for transcripts not associated to any annotated CDS the
  #' position of the P-site with respect to the start and the stop codon is set to
  #' NA. Optionally, additional columns reporting the three nucleotides covered by
  #' the P-site, the A-site and the E-site are attached, based on FASTA files or
  #' BSgenome data packages containing the transcript nucleotide sequences.
  "
}

psite_info <- function(data, offset, site = NULL, fastapath = NULL,
                       fasta_genome = TRUE, refseq_sep = NULL, bsgenome = NULL,
                       gtfpath = NULL, txdb = NULL, dataSource = NA,
                       organism = NA, granges = FALSE) {

  if(!(all(site %in% c("psite", "asite", "esite"))) & length(site) != 0){
    cat("\n")
    stop("parameter site must be either NULL, \"psite\", \"asite\", \"esite\" or a combination of the three strings \n\n")
  } else {
    if(length(site) != 0 & length(fastapath) == 0 & length(bsgenome) == 0){
      cat("\n")
      stop("parameter site is specified but both fastapath and bsgenome are missing \n\n")
    }
  }

  if(length(site) != 0){
    if(((length(fastapath) != 0 & (fasta_genome == TRUE | fasta_genome == T)) |
        length(bsgenome) != 0) &
       length(gtfpath) == 0 & length(txdb) == 0){
      cat("\n")
      stop("genome annotation file not specified (both GTF path and TxDb object are missing)\n\n")
    }

    if(length(fastapath) != 0 & length(bsgenome) != 0){
      cat("\n")
      warning("both fastapath and bsgenome are specified. Only fastapath will be considered\n")
      bsgenome = NULL
    }

    if(length(gtfpath) != 0 & length(txdb) != 0){
      cat("\n")
      warning("both gtfpath and txdb are specified. Only gtfpath will be considered\n")
      txdb = NULL
    }

    if((length(gtfpath) != 0 | length(txdb) != 0) &
       ((length(fastapath) == 0 & length(bsgenome) == 0) |
        (length(fastapath) != 0 & (fasta_genome == FALSE | fasta_genome == F)))){
      cat("\n")
      warning("a genome annotation file is specified but no sequences from genome assembly are provided\n")
    }

    if(length(gtfpath) != 0 | length(txdb) != 0){
      if(length(gtfpath) != 0){
        path_to_gtf <- gtfpath
        txdbanno <- GenomicFeatures::makeTxDbFromGFF(file = path_to_gtf, format = "gtf", dataSource = dataSource, organism = organism)
      } else {
        if(txdb %in% rownames(installed.packages())){
          library(txdb, character.only = TRUE)
        } else {
          source("https://bioconductor.org/biocLite.R")
          biocLite(txdb, suppressUpdates = TRUE)
          library(txdb, character.only = TRUE)
        }
        txdbanno <- get(txdb)
      }
    }

    if(length(fastapath) != 0 | length(bsgenome) != 0){
      if(length(fastapath) != 0) {
        if(fasta_genome == TRUE | fasta_genome == T){
          temp_sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(temp_sequences) <- tstrsplit(names(temp_sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
          exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names = TRUE))
          exon <- as.data.table(exon[unique(names(exon))])
          sub_exon_plus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "+"]
          sub_exon_minus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "-"
                                 ][, new_end := Biostrings::width(temp_sequences[as.character(seqnames)]) - start + 1
                                   ][, new_start := Biostrings::width(temp_sequences[as.character(seqnames)]) - end + 1]

          seq_dt_plus <- sub_exon_plus[, nt_seq := "emp"
                                       ][, nt_seq := as.character(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                                                     start = start,
                                                                                     end = end))
                                         ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]

          revcompl_temp_sequences <- reverseComplement(temp_sequences)
          seq_dt_minus <- sub_exon_minus[, nt_seq := "emp"
                                         ][, nt_seq := as.character(Biostrings::subseq(revcompl_temp_sequences[as.character(seqnames)],
                                                                                       start = new_start,
                                                                                       end = new_end))
                                           ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]

          sequences <- Biostrings::DNAStringSet(c(seq_dt_plus$seq, seq_dt_minus$seq))
          names(sequences) <- c(unique(sub_exon_plus$group_name), unique(sub_exon_minus$group_name))
        } else {
          sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(sequences) <- tstrsplit(names(sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
        }
      } else {
        if(bsgenome %in% installed.genomes()){
          library(bsgenome, character.only = TRUE)
        } else {
          source("http://www.bioconductor.org/biocLite.R")
          biocLite(bsgenome, suppressUpdates = TRUE)
          library(bsgenome, character.only = TRUE)
        }
        sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome), txdbanno, use.names=T)
      }
    }
  }

  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    suboff <- offset[sample == n, .(length,corrected_offset_from_3)]
    cat("1. adding p-site position\n")
    dt[suboff,  on = 'length', psite := i.corrected_offset_from_3]
    dt[, psite := end3 - psite]
    setcolorder(dt,c("transcript", "end5", "psite", "end3", "length", "cds_start", "cds_stop"))
    dt[, psite_from_start := psite - cds_start
       ][cds_stop == 0, psite_from_start := 0]
    dt[, psite_from_stop := psite - cds_stop
       ][cds_stop == 0, psite_from_stop := 0]
    cat("2. adding transcript region\n")
    dt[, psite_region := "5utr"
       ][psite_from_start >= 0 & psite_from_stop <= 0, psite_region := "cds"
         ][psite_from_stop > 0, psite_region := "3utr"
           ][cds_stop == 0, psite_region := NA]
    if(length(site) != 0){
      cat("3. adding nucleotide sequence(s)\n")
      if("psite" %in% site){
        dt[, p_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite,
                                                             end = dt$psite + 2))]
      }
      if("asite" %in% site){
        dt[, a_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite + 3,
                                                             end = dt$psite + 5))]
      }
      if("esite" %in% site){
        dt[, e_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite - 3,
                                                             end = dt$psite - 1))]
      }
    }

    setorder(dt, transcript, end5, end3)

    if (granges == T | granges == TRUE) {
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

    data[[n]] <- dt
  }

  if (granges == T | granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }

  return(data)
}
