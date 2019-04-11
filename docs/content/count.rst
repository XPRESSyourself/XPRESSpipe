############################
Read Quantification
############################

===============================
Quantifying and Collating Reads
===============================
| In order to quantify aligned reads, they must be counts to a reference transcriptome. This will tell you in relative terms how much of each transcript is expressed in a system. The following sub-module will perform this quantification, as well as compile all sample quantifications into a single data matrix for downstream use.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe count --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of SAM files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-r \<path\>, --reference \<path\>`
     - Path to parent organism reference directory
   * - :data:`-g \</path/transcript.gtf\>, --gtf \</path/transcript.gtf\>`
     - Path and file name to GTF used for alignment quantification

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-e <experiment_name>`, :data:`--experiment <experiment_name>`
     - Experiment name
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Count ribosome profiling alignments:**
| - Input points to directory with SAM alignment files that are sorted by name
| - An experiment name is provided to name the final data matrix
| - Reads are quantified only to coding genes and are not counted if mapping to the first x nucleotides of each transcript exon 1 (x being the value provided for truncation when initially creating the reference files)

.. code-block:: shell

  $ xpresspipe count -i riboprof_out/alignments/ -o riboprof_out/ -r se_reference/ -g se_reference/transcripts_codingOnly_truncated.gtf -e se_test

| **Example 2 -- Count paired-end alignments:**
| - Input points to directory with SAM alignment files that are sorted by name
| - An experiment name is not provided and a default name is given to the data matrix using datatime
| - Reads are quantified to the entire transcriptome (coding and non-coding, no truncation)

.. code-block:: shell

  $ xpresspipe count -i pe_out/alignments/ -o pe_out/ -r pe_reference/

----------------------
Future Implementations
----------------------
| Add arguments for -s strandedness so other kits are usable with pipeline (default: False for TCGA)
| Add arguments for -m intersection (d)efault: intersection-nonempty for TCGA)
