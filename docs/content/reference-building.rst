###################
Curating References
###################
| In order to quantify transcription on a transcript to transcript basis, individual reads called during sequencing must be mapped to the genome. While there are multiple alignment software packages available, XPRESSpipe uses a current version of `STAR <https://github.com/alexdobin/STAR>`_ to perform this step in transcription quantification for several reasons:
| - Performance: While computationally greedy (a human genome alignment requires upwards of 30 Gb RAM), the `performance and accuracy is superior to the majority of other splice aware aligners currently available <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/>`_
| - Splice Junction Aware: STAR is capable of mapping reads spanning a splice junction, where more traditional packages, such as Bowtie, are incapable of doing so and are better suited for tasks such as genome alignment.
| - Standard: The foundation of the pipeline used in XPRESSpipe is based in the `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ standards for RNAseq alignment. This method utilizes a 2-pass alignment program, where alignment is performed to identify splice junctions, these junctions are considered in the reference genome, and reads are re-aligned with this new reference. While more time-intensive, the practice provides the user with a more thorough quantification of their RNA sequencing.

=================================
XPRESSpipe Reference Requirements
=================================
| An XPRESSpipe compatible reference directory must meet some requirements:
| - All chromosomal genome fasta files are in the parent reference directory
| - A sub-directory, named :data:`genome`, contains the STAR reference files. If :data:`createReference` is used to curate the reference, and the parent reference directory was provided as output location, this formatted will be handled by XPRESSpipe
| - A transcript reference (GTF), is located in the reference parent directory and is named :data:`transcripts.gtf`. If coding-only or truncated reference GTFs are used, they should also be in this directory (:data:`truncate` will handle file naming and formatting so long as the output location is specified as this parent directory)
| - If :data:`metagene` analysis is being performed, properly flattened transcript references are located in the parent reference directory. If :data:`makeFlat` was used with this parent reference directory as output location, proper naming and formatted will be handled by XPRESSpipe
| **A completed reference directory can be created that follows these requirements by creating a directory, placing the transcripts.gtf and chromosome genomic fasta files in the parent directory and running :data:`curateReference` as described below

==========================
STAR Reference Curation
==========================
| The following creates a STAR reference compatible with XPRESSpipe. These files are output in a directory called :data:`genome` in the specified :data:`--output` directory.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe createReference --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-o \<str\>, --output \<str\>`
     - Path to output directory
   * - :data:`-f \<str\>, --fasta \<str\>`
     - Path to genome fasta files (file names should end in .fa, .fasta, or .txt and no other files should exist in the directory with similar extensions)
   * - :data:`-g \<str\>, --gtf \<str\>`
     - Path and file name to transcript reference file names 'transcripts.gtf'

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-t <int>, --threads <int>`
     - Specify number of threads to use (default: :data:`8`)
   * - :data:`--sjdbOverhang \<int\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases

-----------
Examples
-----------
| **Example 1 -- Create a single-end sequencing reference:**
| - Paths to output and location of genome fasta files for each chromosome are provided, as well as path and file name to transcripts.gtf file
| - Default number of threads are used for preparing reference

.. code-block:: shell

  $ xpresspipe createReference -o /path/to/reference/ -f /path/to/reference/ -g /path/to/reference/transcripts.gtf -sjdbOverhang 49

| **Example 2 -- Create a paired-end sequencing reference:**
| - 12 threads are specified for reference creation
| - The as 2x100bp paired-end sequencing was used, the default value for :data:`--sjdbOverhang` of :data:`100` is appropriate in this case

.. code-block:: shell

  $ xpresspipe createReference -o /path/to/reference/ -f /path/to/reference/ -g /path/to/reference/transcripts.gtf -t 12

============================================
Transcript Reference Curation and Truncation
============================================
| At times, quantification of transcripts to coding-only transcripts or to a modified transcript reference is desirable. Below are two examples:
| 1. As ribosomal RNA (rRNA) contamination is common in RNAseq, even when a depletion step was performed prior to library preparation, it is sometimes desirable to not count these and other non-coding RNAs in the quantification and analysis.
| 2. During ribosome profiling library preparation, a 5' transcript bias is common, regardless of library preparation method. It has therefore been suggested to `exclude the first 45-50 nucleotides of each transcript from quantification <https://www.cell.com/cms/10.1016/j.celrep.2016.01.043/attachment/257faf34-ff8f-4071-a642-bfdb531c75b8/mmc1>`_

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe truncate --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-g \<str\>, --gtf \<str\>`
     - Path and file name to transcript reference file (GTF) to process

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-t <int>, --truncate_amount <int>`
     -  Number of nucleotides to truncate from the 5' end of each transcript (default: :data:`45`)
   * - :data:`-c, --create_refFlats`
     - Provide flag to output refFlat files for each transcript reference created

-----------
Examples
-----------
| **Example 1 -- Create coding-only and truncated references:**
| - Creates a coding only GTF reference file and a truncated coding-only reference file
| - Truncates the first 50 nucleotides from the first exon of every transcript

.. code-block:: shell

  $ xpresspipe truncate -g /path/to/reference/transcripts.gtf -t 50 -c

============================================
Flatten Transcript References
============================================
| Certain analysis platforms require a RefFlat transcript prediction file. These files can be created with the following command.
| XPRESSpipe uses the `UCSC-GTFtoGenePred <https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html>`_ package to perform these conversions.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe makeFlat --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<str\>, --input \<str\>`
     - Path where input transcripts*.gtf files are found

-----------
Examples
-----------
| **Example 1 -- Create refFlat files:**
| - Creates a refFlat-formatted file for each GTF file in the given input directory

.. code-block:: shell

  $ xpresspipe makeFlat -i /path/to/reference/

============================================
Perform Full Reference Curation
============================================
| The following will create a XPRESSpipe-formatted reference directory containing all STAR reference files and transcript references needs for quantification and meta-analysis.
| A parent reference directory containing the transcripts.gtf file and all chromosomal genome fasta files must be present

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe curateReference --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-o \<str\>, --output \<str\>`
     - Path to output directory
   * - :data:`-f \<str\>, --fasta \<str\>`
     - Path to genome fasta files (file names should end in .fa, .fasta, or .txt and no other files should exist in the directory with similar extensions)
   * - :data:`-g \<str\>, --gtf \<str\>`
     - Path and file name to transcript reference file names 'transcripts.gtf'

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-t <int>, --threads <int>`
     - Specify number of threads to use (default: :data:`8`)
   * - :data:`--sjdbOverhang \<int\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`-t <int>, --truncate_amount <int>`
     -  Number of nucleotides to truncate from the 5' end of each transcript (default: :data:`45`)

-----------
Examples
-----------
| **Example 1 -- Create XPRESSpipe-formatted reference for single-end alignment:**
| - Creates a star reference for single-end read mapping (1x50bp reads)
| - Outputs coding-only and truncated coding-only transcripts reference GTFs
| - Truncates the first 50 nucleotides from the first exon of every transcript
| - Creates a refFlat-formatted file for each GTF file in the given input directory

.. code-block:: shell

  $ xpresspipe curateReference -o $SE_REF/ -f $SE_REF/ -g $SE_REF/transcripts.gtf -t 50 -m 10 --sjdbOverhang 49

| **Example 2 -- Create refFlat files:**
| - Creates a star reference for paired-end read mapping (2x100bp reads)
| - Outputs coding-only and truncated coding-only transcripts reference GTFs
| - Creates a refFlat-formatted file for each GTF file in the given input directory

.. code-block:: shell

  $ xpresspipe curateReference -o $PE_REF/ -f $PE_REF/ -g $PE_REF/transcripts.gtf -m 10
