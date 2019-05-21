# Read quantification
[X] Add count_cufflinks
[X] Add cufflinks aggregator for table
[X] No normalization needed
[ ] Add masking creation in curation
  - cat Homo_sapiens.GRCh38.95.gtf | awk '$1 == "MT"' | ct
  - cat Homo_sapiens.GRCh38.95.gtf | grep -i "rrna\|trna\|misc_RNA" | ct
  - Add masking to htseq?
[X] Add options to arguments
[X] Edit main for needed changes
[ ] Check for macOS and quartz on install, if no, print error
[ ] Add cufflinks info to docs
[X] Pipeline test

# Truncator
[ ] Find exon space per transcript and anything less than or equal to \_5prime + \_3prime remove all associated gene records


# Analysis

## Metagene


## Periodicity



## Complexity



## Read distribution



# Pipeline

## SE



## PE  



## riboprof
