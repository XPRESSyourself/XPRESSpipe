###########
Trim
###########

| **Trimming is a necessary part of RNAseq data processing due to the technological limitations described below:**
| - Inherent in RNAseq library creation, RNA is fragmented and adaptor sequences are ligated to the sequence. These adaptors include information such as sample batch and act as a primer for the sequencer to recognize the fragment as something to analyze. However, these adaptors, once sequenced, prevent alignment to a reference as large chunks of the fragment are synthetic sequence not found in each transcript.
| - A sequencer's job is to read a fragment base by base and determine the nucleotide species each step of the way. While the technology has greatly improved over the years, a probability of error remains. Miscalled bases can prevent proper alignment of the sequenced fragment to the reference. Therefore, it is important for low confidence base calls to be trimmed from each read.
|
| Trimming is performed by `fastp <https://github.com/OpenGene/fastp>`_

================
Arguments
================
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe trim --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i INPUT, --input INPUT`
     - Path to input directory -- if paired-end, file names should be exactly the same except for :data:`r1/r2.fastq` or similar suffix
   * - :data:`-o OUTPUT, --output OUTPUT`
     - Path to output directory

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-a \<list\> [\<list\> ...]`, :data:`--adaptor \<list\> [\<list\> ...]`
     - Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if :data:`None` are provided, software will attempt to auto-detect adaptors -- if "POLYX" is provided as a single string in the list, polyX adaptors will be trimmed. If you want to auto-detect adaptors in for paired-end reads, provide :data:`None` twice
   * - :data:`-q \<int\>, --quality \<int\>`
     - PHRED read quality threshold (default: :data:`28`)
   * - :data:`--min_length \<int\>`
     - Minimum read length threshold to keep for reads (default: :data:`18`)
   * - :data:`-m <int>, --max_processors <int>`
     - Number of max processors to use for tasks (default: Max)

================
Examples
================
| **Example 1 -- Trim ribosome profiling (or single-end) sequence data using default preferences:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i riboprof_test/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o riboprof_out/`
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/

| **Example 2 -- Predict adaptor and trim ribosome profiling (or single-end) sequence data:**
| - A minimum read length of 22 nucleotides after trimming is required in order to keep the read
| - A maximum or 6 processors can be used for the task
| - The :data:`--adaptors` argument was not passed, so an attempt to discover adaptor sequences will be made (this is not always the most efficient or thorough method of trimming and providing the adaptor sequences is recommended)

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6

| **Example 3 -- Pass explicit adaptor trim ribosome profiling (or single-end) sequence data:**
| - The default minimum read length threshold will be used
| - The maximum number of processors will be used by default
| - The :data:`--adaptors` argument was passed, so adaptor sequences will trimmed explicitly

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT

| **Example 4 -- Predict adaptor and trim paired-end sequence data:**
| - The :data:`--adaptors` argument was passed as :data:`None None`, so an attempt to discover adaptor sequences will be made for paired-end reads. The :data:`-a None None` syntax is essential for :data:`trim` to recognize the reads as paired-end

.. code-block:: shell

  $ xpresspipe trim -i pe_test/ -o pe_out/ -a None None

| **Example 5 -- Pass explicit adaptor and trim paired-end sequence data:**
| - The :data:`--adaptors` argument was passed, so adaptor sequences will trimmed explicitly

.. code-block:: shell

  $ xpresspipe trim -i pe_test/ -o pe_out/ -a ACACTCTTTCCCTACACGACGCTCTTCCGATC GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

| **Example 6 -- Trim single-end sequence data of polyX adaptors:**
| - The :data:`--adaptors POLYX` argument was passed, so adaptor sequences will trimmed of polyX sequences

.. code-block:: shell

  $ xpresspipe trim -i se_test/ -o se_out/ -a POLYX
