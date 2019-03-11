############################
Normalize
############################

| Note: Sample and batch normalization can be performed in a single command. If this is done, batch normalization will be performed following sample normalization.

================
Sample Normalize
================
| Due to inherent biases in RNAseq samples (most commonly, different amounts of total RNA per sample in a given lane), samples must be normalized to obtain an accurate representation of transcription per sample. Additional normalization can be performed to normalize for transcript length ("per kilobase million") as longer transcripts will naturally create more fragments mapping to a given gene, thus potentially making 1 transcript appear as many when quantified.
|
| The following equations summarize different way to normalize samples for RNAseq:
| Perform reads per million sample normalization on RNAseq data
| :math:`RPM = \frac{\#\ number\ reads\ per\ gene\ x\ 1e6}{\#\ mapped\ reads\ per\ sample}`
| Perform reads per kilobase million sample normalization on RNAseq data
| :math:`RPKM = \frac{\#\ number\ reads\ per\ gene\ x\ 1e6\ x\ 1e3}{\#\ mapped\ reads\ per\ sample\ x\ gene\ length\ (bp)}`
| Perform fragments per kilobase million sample normalization on RNAseq data
| :math:`FPKM = \frac{\#\ number\ fragments\ per\ gene\ x\ 1e6\ x\ 1e3}{\#\ mapped\ fragments\ per\ sample\ x\ gene\ length\ (bp)}`
|
| Translation efficiency normalization can be performed within XPRESStools
|
| Assumptions:
|   - R is installed on your machine and is in your $PATH
|   - All input files are tab-delimited (with .txt or .tsv suffix)

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe normalizeMatrix --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-d \<path/filename.tsv\>, --data \<path/filename.tsv\>`
     - Path and file name of expression counts matrix

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`--method \<RPM, RPKM, FPKM, LOG\>`
    - Normalization method to perform (options: "RPM", "RPKM", "FPKM", "LOG") -- if using either RPKM or FPKM, a GTF reference file must be included
  * - :data:`-g \</path/transcripts.gtf\>, --gtf \</path/transcripts.gtf\>`
    - Path and file name to reference GTF
  * - :data:`--batch \</path/filename.tsv\>`
    - Include path and filename of dataframe with batch normalization parameters

-----------
Examples
-----------
| **Example 1 -- Perform RPKM normalization on single-end RNAseq data:**

.. code-block:: shell

  $ xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.tsv --method RPKM -g se_reference/transcripts_coding_truncated.gtf

=====================
Batch Normalize
=====================
| When multiple people perform library preparation, or when libraries are prepared on different days, this can lead to inherent biases in count distributions between batches of samples. It is therefore necessary to normalize these effects when appropriate.
|
| Assumptions:
|   - R is installed on your machine and is in your $PATH
|   - All input files are tab-delimited (with .txt or .tsv suffix)

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe normalizeMatrix --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-d \<path/filename.tsv\>, --data \<path/filename.tsv\>`
     - Path and file name of expression counts matrix

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`--method \<RPM, RPKM, FPKM, LOG\>`
    - Normalization method to perform (options: "RPM", "RPKM", "FPKM", "LOG") -- if using either RPKM or FPKM, a GTF reference file must be included
  * - :data:`-g \</path/transcripts.gtf\>, --gtf \</path/transcripts.gtf\>`
    - Path and file name to reference GTF
  * - :data:`--batch \</path/filename.tsv\>`
    - Include path and filename of dataframe with batch normalization parameters

-----------
Examples
-----------
| **Example 1 -- Perform batch normalization on RNAseq data:**

.. ident with TABs
.. code-block:: python

  > batch = pd.read_csv('./riboprof_out/counts/batch_info.tsv', sep='\t', index_col=0)
  > batch
    Sample  Batch
  0 s1      batch1
  1 s2      batch2
  2 s3      batch1
  3 s4      batch2

.. code-block:: shell

  $ xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.tsv --batch riboprof_out/counts/batch_info.tsv
