############################
Ribosome Profiling Pipeline
############################
| The following pipeline will pre-process, align, and quality check ribosome profiling samples using the sub-modules discussed in earlier chapters. For more detailed information concerning these steps, please refer to the appropriate chapter.

| Steps:
.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Function
     - Description
   * - :data:`run_trim()`
     - Trims reads of adaptors, low-quality base calls, and reads that are too short
   * - :data:`run_seRNAseq()`
     - Aligns each samples reads using a two-pass splice-aware alignment; formats output files for downstream use
   * - :data:`create_bed()` -- optional
     - Outputs bed-formatted file
   * - :data:`create_bigwig()` -- optional
     - Outputs bigwig-formatted file
   * - :data:`count_reads()`
     - Quantitates aligned reads per each sample
   * - :data:`collect_counts()`
     - Creates a table of counts for each sample
   * - :data:`run_normalization()`
     - Normalizes count table based on user-defined input
   * - :data:`make_readDistributions()`
     - Runs fastqc on each sample and generates a summary for read distributions of each sample
   * - :data:`make_periodicity()`
     - Generates a summary periodicity report for each ribosome profiling footprint sample
   * - :data:`make_metagene()`
     - Generates a summary meta-gene report for each sample
   * - :data:`get_summary()`
     - Runs multiqc to generate summary reports for steps in the pipeline

-----------
Examples
-----------
| **Example 1 -- Run pipeline on ribosome profiling sample files**

.. ident with TABs
.. code-block:: python

  $ xpresspipe riboprof -i riboprof_test -o ribopipe_out -r se_reference -t CODING_TRUNCATED -e riboprof_test -a CTGTAGGCACCATCAAT --method RPKM --sjdbOverhang 49
