############################
Alignment
############################
| In order to quantify transcription on a transcript to transcript basis, individual reads called during sequencing must be mapped to the genome. While there are multiple alignment software packages available, XPRESSpipe uses a current version of `STAR <https://github.com/alexdobin/STAR>`_ to perform this step in transcription quantification for several reasons:
| - Performance: While computationally greedy (a human genome alignment requires upwards of 30 Gb RAM), the `performance and accuracy is superior to the majority of other splice aware aligners currently available <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/>`_
| - Splice Junction Aware: STAR is capable of mapping reads spanning a splice junction, where more traditional packages, such as Bowtie, are incapable of doing so and are better suited for tasks such as genome alignment.
| - Standard: The foundation of the pipeline used in XPRESSpipe is based in the `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ standards for RNAseq alignment. This method utilizes a 2-pass alignment program, where alignment is performed to identify splice junctions, these junctions are considered in the reference genome, and reads are re-aligned with this new reference. While more time-intensive, the practice provides the user with a more thorough quantification of their RNA sequencing.

============================
Single-End RNAseq Alignment
============================
| The following runs single-end reads alignment using the specified XPRESSpipe-formatted reference directory.
| Notes:
| - For the :data:`--sjdbOverhang` argument, the same value should be entered that was used when creating the STAR reference files.
| - Ribosome profiling data can be run as a single-end RNAseq

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe align --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<str\>, --input \<str\>`
     - Path to input directory
   * - :data:`-o \<str\>, --output \<str\>`
     - Path to output directory
   * - :data:`-t \<str\>, --type \<str\>`
     - Sequencing type ("SE" for single-end, "PE" for paired-end)
   * - :data:`-r \<str\>, --reference \<str\>`
     - Path to parent organism reference directory

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`--output_bigwig`
     - Include flag to output bigwig files for each aligned file
   * - :data:`-m <int>, --max_processors <int>`
     - Number of max processors to use for tasks (default: Max)
   * - :data:`--sjdbOverhang \<int\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases

-----------
Examples
-----------
| **Example 1 -- Single-end RNAseq alignment:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i /path/to/input/files/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o riboprof_out/`
| - :data:`--type` is specified as 'SE' and path to parent reference directory is provided
| - The value for :data:`--sjdbOverhang` used in reference creation is provided. Failure to do so will trigger an error
| - BED and BIGWIG files will be output in their own directories in :data:`output`
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe align -i /path/to/input/files/ -o riboprof_out/ -t SE -r /path/to/reference/ --sjdbOverhang 49 --output_bed --output_bigwig

============================
Paired-End RNAseq Alignment
============================
| The following runs paired-end reads alignment using the specified XPRESSpipe-formatted reference directory.
| Notes:
| - For the :data:`--sjdbOverhang` argument, the same value should be entered that was used when creating the STAR reference files.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe align --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<str\>, --input \<str\>`
     - Path to input directory
   * - :data:`-o \<str\>, --output \<str\>`
     - Path to output directory
   * - :data:`-t \<str\>, --type \<str\>`
     - Sequencing type ("SE" for single-end, "PE" for paired-end)
   * - :data:`-r \<str\>, --reference \<str\>`
     - Path to parent organism reference directory

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`--output_bigwig`
     - Include flag to output bigwig files for each aligned file
   * - :data:`-m <int>, --max_processors <int>`
     - Number of max processors to use for tasks (default: Max)
   * - :data:`--sjdbOverhang \<int\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases

-----------
Examples
-----------
| **Example 1 -- Paired-end RNAseq alignment:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i pe_test/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o pe_out/`
| - :data:`--type` is specified as 'PE' and path to parent reference directory is provided
| - The value for :data:`--sjdbOverhang` used in reference creation is provided. Failure to do so will trigger an error. In this case, since the reference was created using default values, the optional flag is not used
| - BED and BIGWIG files are not output
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe align -i /path/to/input/files/ -o riboprof_out -t PE -r /path/to/reference/
