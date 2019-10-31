###############
Updates
###############

============
v0.3.1
============
- Fix BAM file threshold for metagene and geneCoverage to avoid OOM errors 
- Turn off BAM file threshold for counting (low memory footprint, so can use all cores available)
- Import openssl library manually in Rperiodicity -- occasionally had trouble finding the library on its own and would error

============
v0.3.0
============
| - Transfers R dependency installs to Anaconda environment load
| - Modified fastq and bam memory factor to optimize resources
| - Rebuilt read distribution module with JuliaLang for super memory efficiency during parallelization
| - Fixed issue where one| -exon genes would not display feature annotations in `geneCoverage` modules
| - Made matplotlib backend calls flexible for HPC usage
| - Made directory checks more thorough
| - Fixed a potential off| -by| -one issue with GTF truncator
| - Updated appropriate tests
| - Updates to documentation
| - Added code of conduct and contributions information

===========
v0.2.4-beta
===========
| - Manuscript submission version
| - Fixed issues with using polyX adaptors
| - Allowed more multi-threading during post-processing of aligned reads to use resources more efficiently
| - Added integrated pipeline tests for Travis CI build to assess pipeline integrity each push
| - Updated install walkthrough video
