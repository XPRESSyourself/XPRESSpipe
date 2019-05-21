# Read quantification
[X] Add count_cufflinks
[X] Add cufflinks aggregator for table
[X] No normalization needed
[X] Add options to arguments
[X] Edit main for needed changes
[X] Check for macOS and quartz on install, if no, print error
[X] Pipeline test
[ ] Add cufflinks info to docs
[ ] Test on HPC


# Alignment
[X] Curate masking index for STAR
  - curate
    - ask if mask true and point to mask file
  - make
    - require path to mask file
[X] Implement masking step in alignment in ALIGN
[X] Test
[ ] Add info to docs
  - ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/ncrna/
[ ] Test on HPC


# Truncator
[ ] Find exon space per transcript and anything less than or equal to \_5prime + \_3prime remove all associated gene records


# Analysis
[ ] Check logic
[ ] Check function

## Metagene


## Periodicity


## Complexity


## Read distribution
[ ] Build own


# Pipeline
## SE
[ ] Run without errors and expected output (parse from log file?)


## PE  
[ ] Run without errors and expected output (parse from log file?)


## riboprof
[ ] Run without errors and expected output (parse from log file?)
