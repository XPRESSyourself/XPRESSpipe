############################
Ribosome Profiling Pipeline
############################
| The following pipeline will pre-process, align, and quality check ribosome profiling samples using the sub-modules discussed in earlier chapters. For more detailed information concerning these steps, please refer to the appropriate chapter.

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory -- if paired-end, file names should be exactly the same except for :data:`r1/r2.fastq` or similar suffix
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-r \<path\>, --reference \<path\>`
     - Path to parent organism reference directory
   * - :data:`-t \<SE or PE\>, --type \<SE or PE\>`
     - Sequencing type ("SE" for single-end, "PE" for paired-end)
   * - :data:`-e <experiment_name>`, :data:`--experiment <experiment_name>`
     - Experiment name

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--two-pass`
     - Use a two-pass STAR alignment for novel splice junction discovery
   * - :data:`-a \<adaptor1 ...\> [\<adaptor1 ...\> ...]`, :data:`--adaptor \<adaptor1 ...\> [\<adaptor1 ...\> ...]`
     - Specify adaptor(s) in list of strings -- for single-end, only provide one adaptor -- if :data:`None` are provided, software will attempt to auto-detect adaptors -- if "POLYX" is provided as a single string in the list, polyX adaptors will be trimmed. If you want to auto-detect adaptors in for paired-end reads, provide :data:`None` twice
   * - :data:`-q \<PHRED_value\>, --quality \<PHRED_value\>`
     - PHRED read quality threshold (default: :data:`28`)
   * - :data:`--min_length \<length_value\>`
     - Minimum read length threshold to keep for reads (default: :data:`18`)
   * - :data:`--deduplicate`
     - Include flag to quantify reads with de-duplication (will search for files with suffix :data:`_dedupRemoved.bam`)
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`-c <method>`, :data:`--quantification_method <method>`
     - Specify quantification method (default: cufflinks; other option: htseq. If using Cufflinks, no downstream sample normalization is required)
   * - :data:`--method \<RPM, RPKM, FPKM, LOG\>`
     - Normalization method to perform (options: "RPM", "TPM", "RPKM", "FPKM") -- if using either TPM, RPKM, or FPKM, a GTF reference file must be included
   * - :data:`--batch \</path/filename.tsv\>`
     - Include path and filename of dataframe with batch normalization parameters
   * - :data:`--sjdbOverhang \<sjdbOverhang_amount\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`--mismatchRatio \<mismatchRatio\>`
     - Alignment ratio of mismatches to mapped length is less than this value. See STAR documentation for more information on setting this parameter
   * - :data:`--seedSearchStartLmax \<seedSearchStartLmax\>`
     - Adjusting this parameter by providing a lower number will improve mapping sensitivity (recommended value = 15 for reads ~ 25 nts). See STAR documentation for more information on setting this parameter
   * - :data:`genome_size`
     - Only needs to be changed if this argument was provided curing reference building AND using a two-pass alignment
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

| Run the following for more details:

.. ident with TABs
.. code-block:: python

  $ xpresspipe riboseq --help

-----------
Examples
-----------
| **Example 1 -- Run pipeline on ribosome profiling sample files**

.. ident with TABs
.. code-block:: python

  $ xpresspipe riboseq \
                -i riboprof_test \
                -o ribopipe_out \
                -r se_reference \
                --gtf se_reference/transcript_LCT.gtf \
                -e riboprof_test \
                -a CTGTAGGCACCATCAAT \
                --method RPKM \
                --sjdbOverhang 49
