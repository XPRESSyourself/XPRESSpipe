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
sub_license <- function() {
  "
  Portions of this module are derived from riboWaltz:

  MIT License

  Copyright (c) [2017] [LAboratory of Translational Architectomics]

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the 'Software'), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
  "
  }

# Build index from GTF for quality control analyses

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org")}

install.packages("stringi", repos = "http://cran.us.r-project.org")
library(stringi)

if ("GenomicFeatures" %in% rownames(installed.packages()) == FALSE) {
  print("Installing GenomicFeatures...")
  BiocManager::install("GenomicFeatures", dependencies=c("Depends", "Imports", "LinkingTo"))
} else {
  print("GenomicFeatures package already installed")
}

library(GenomicFeatures)

if ("data.table" %in% rownames(installed.packages()) == FALSE) {
  print("Installing data.table...")
  install.packages("data.table", repos = "http://cran.us.r-project.org")
} else {
  print("data.table package already installed")
}
library(data.table)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
GTF <- args[1]
OUTPUT <- args[2]

# Create flat GTF file for quality control modules with transcript id, gene name, and feature lengths
flatten_gtf <- function(
  gtf_file) {

  if (endsWith(gtf_file, '.gtf')) {

    # Read in file to Genomic Alignments data.frame and report number of records
    print(paste('Reading ', gtf_file, '...', sep=''))
    gtf <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file, format = 'gtf')

    # Get features
    gene <- suppressWarnings(GenomicFeatures::transcriptsBy(gtf, "gene"))
    exon <- suppressWarnings(GenomicFeatures::exonsBy(gtf, by = "tx",use.names=T))
    utr5<- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(gtf,use.names=T))
    cds <- suppressWarnings(GenomicFeatures::cdsBy(gtf, by = "tx", use.names=T))
    utr3<- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(gtf,use.names=T))

    # Convert features to dataframes
    gene <- as.data.table(gene)
    exon <- as.data.table(exon[unique(names(exon))])
    utr5 <- as.data.table(utr5[unique(names(utr5))])
    cds <- as.data.table(cds[unique(names(cds))])
    utr3 <-as.data.table(utr3[unique(names(utr3))])

    # Parse out relevant information, with transcript_id being the key
    gene_name <- gene[, c('tx_name','group_name')]
    names(gene_name) <- c("transcript", "gene")
    l_transcript <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
    l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
    l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
    l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]

    # Merge records
    print('Collecting flat reference information...')
    merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
    flat_gtf  <-  Reduce(merge_allx, list(gene_name, l_transcript, l_utr5, l_cds, l_utr3))
    flat_gtf[is.na(flat_gtf)] <- 0

    rows <- dim(flat_gtf)[1]
    print(paste(rows, 'transcript records compiled...', sep=' '))

    return(flat_gtf)

  } else {

    print(paste(gtf_file, 'is an invalid file type. File must be a GTF file with the .gtf suffix', sep=' '))

  }

  }

# MAIN
reference <- flatten_gtf(GTF)
write.table(as.data.frame(reference), file=OUTPUT, sep='\t', col.names=NA)
