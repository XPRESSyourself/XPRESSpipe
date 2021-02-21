###############
Updates
###############

========================
v0.6.2
========================
| - Added :data:`--suppress_version_check` flag to enable use of XPRESSpipe without internet access
| - Added :data:`--smoothen` flag to any module that uses the :data:`geneCoverage` sub-module. By default, a sliding window will not be used to smoothen the geneCoverage plots. If provided, a rolling window set at 20 will be used to smoothen the plots.

========================
v0.6.1
========================
| - Add flag during curation steps to allow of UCSC/refseq GTFs during GTF modification steps (truncation, etc.)
|   Usage: Provide the :data:`--ucsc_format` flag to the :data:`curateReference` or :data:`modifyGTF` sub-modules. These modifications in format only apply to XPRESSpipe GTF truncation features. Any formatting errors with the GTF file that pertain to alignment, counting, etc. dependencies will need to be addressed by the user.
| - Fixed error in XPRESSpipe interface with XPRESSplot's convert_names function where XPRESSpipe did not read in first column of table as index

========================
v0.6.0
========================
| - Minor modification to instructions on how to install XPRESSpipe and use its conda environment on a supercomputing node.
| - Removed version specifications for conda environment setup to ease install (fixes issues in a better way than the solution from v0.5.0)
| - Moved required riboWaltz functions to XPRESSpipe as installation has been recurrently problematic

========================
v0.5.0
========================
| - Fixed issue where genome size calculation would round up and miscalculate genome_size parameter for STAR.
| - Added fastp_lite for removal of 3' internal UMIs (generally takes ~1 min per RNA-seq sample with about 30 million reads)
|   - For example:
|     5'-read-spacer-UMI-adapter-3'
|   - Integrating this addition into options and trimming
| - Updated MANIFEST file to reliably copy R and Julia scripts to executing folder for XPRESSpipe
| - Updated command builder (:data:`xpresspipe build`) to include recent additions
| - Updated requirements to prevent issue where solved environment required to install :data:`R 3.5.1` or greater would create an error where :data:`samtools markdup` would freeze
| - Frequently, :data:`R 3.5.1` or greater would create library linking error to :data:`stringi`, causing :data:`GenomicFeatures` to not function. Added to :data:`RbuildIndex.r` to reinstall :data:`stringi`, which appears to clear up the issue.


========================
v0.4.4
========================
| - Fixed issue with string catenation during UMI fastp call where UMI length was not properly forced to a string

========================
v0.4.3
========================
| - Fixed issue with `metagene` where parallelization overloaded memory and resulted in OOM errors. Fixed by making memory thresholding slightly more strict.

========================
v0.4.2
========================
| - Fixed `convert_names` xpressplot call

========================
v0.4.1
========================
| - Introduced some restrictions to dependency versions. Some newer versions of dependencies were acting problematic. Will try to figure out how to allow for current versions of these dependencies to be used
| - Fixed plotting issue with periodicity plots

========================
v0.4.0
========================
| - Introduced rRNA depletion during alignment step (previously could only do so during the quantification step)
| - Expanded periodicity analysis to cover more holistic P-site analysis (report codon usage stats). The submodule previously called :data:`periodicity` is now called by :data:`p_sites`
| - Allow for setting upper limit threshold for read length during pre-processing reads and P-site analysis (previously only a lower limit was available)
| - All documentation associated with this changes has been updated.

============
v0.3.1
============
| - Fix BAM file threshold for metagene and geneCoverage to avoid OOM errors
| - Turn off BAM file threshold for counting (low memory footprint, so can use all cores available)
| - Import openssl library manually in Rperiodicity -- occasionally had trouble finding the library on its own and would error

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
