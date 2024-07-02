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
  Annotation data table.

  This function generates transcript basic annotation data tables starting from
  GTF files or TxDb objects. Annotation data tables include a column named
  transcript reporting the name of the reference transcripts and four
  columns named l_tr, l_utr5, l_cds and l_utr3
  reporting the length of the transcripts and of their annotated 5' UTRs, CDSs
  and 3' UTRs, respectively. Please note: if a transcript region is not
  annotated its length is set to 0.
  "
}

library(data.table)
library(GenomicFeatures)
#library(txdbmaker, character.only = TRUE)

create_annotation <-  function(gtfpath = NULL, txdb = NULL, dataSource = NA, organism = NA) {

  if(length(gtfpath) != 0 & length(txdb) != 0){
    cat("\n")
    warning("gtfpath and txdb are both specified. Only gtfpath will be considered\n")
    txdb = NULL
  }

  if(length(gtfpath) == 0 & length(txdb) == 0){
    cat("\n")
    stop("neither gtfpath nor txdb is specified\n\n")
  }

  if(length(gtfpath) != 0){
    path_to_gtf <- gtfpath
    print(path_to_gtf)
    txdbanno <- txdbmaker::makeTxDbFromGFF(file = path_to_gtf, format = "gtf", dataSource = dataSource, organism = organism)
  } else {
    txdbanno <- get(txdb)
  }

  exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx",use.names=T))
  utr5<- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdbanno,use.names=T))
  cds <- suppressWarnings(GenomicFeatures::cdsBy(txdbanno, by = "tx", use.names=T))
  utr3<- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdbanno,use.names=T))
  exon <- as.data.table(exon[unique(names(exon))])
  utr5 <- as.data.table(utr5[unique(names(utr5))])
  cds <- as.data.table(cds[unique(names(cds))])
  utr3 <-as.data.table(utr3[unique(names(utr3))])

  anno_df <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
  l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
  l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
  l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]

  merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
  anno_df  <-  Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0

  return(anno_df)
}
