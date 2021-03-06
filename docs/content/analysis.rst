.. _analysis_link:

Analysis
##############################

=================================
Differential Expression Analysis
=================================
| Differential Expression analysis allows one to determine significantly enriched or depleted genes between two conditions. XPESSpipe acts as a wrapper for `DESeq2 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/>`_. Please refer to its `documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_ for more information.
| NOTE: If intending to use the :data:`diffxpress` sub-module, you need to have used :data:`--quantification_method htseq` during read quantification as `DESeq2 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/>`_ requires integer count data.

|
| Requirements:
|   - R is installed on your machine and is in your $PATH (this should be handled in the installation)
|   - All input files are tab-delimited (with .txt or .tsv suffix)
|   - Design formula does not include the tilde (~) and there are no spaces


---------------------
Sample Factor Files
---------------------
| Different factors to be evaluated in the differential expression analysis should each be denoted as a separate factor column in the :data:`sample_info` file. For example, if you were evaluating a experimental vs control experiment for RNA-sequencing, you would provide a :data:`sample_info` file as follows:

.. ident with TABs
.. code-block:: python

  **sample_info.txt**

    Sample      Condition
    s1_rna      a_WT
    s2_rna      a_WT
    s3_rna      b_EXP
    s4_rna      b_EXP

| Your base (denominator) parameter in a given factor column in the :data:`sample_info` file must be first alphabetically. In the case provided above, we want to compare the experimental condition *VS* the wild-type control condition, however these labels are not alphabetical. In this case, you can append letters to the beginning to force alphabetical order. For example, if you performed a :data:`experiment` vs :data:`wild-type` experiment, you would need to use the labels :data:`b_experiment` vs :data:`a_wild-type` to force a :data:`b_experiment` / :data:`a_wild-type` comparison.

| If we want to consider additional factors, such as translation efficiency of  footprint vs RNA-sequence samples for ribosome profiling, these should be included as additional factor columns in the :data:`sample_info` file. Since we want to perform another comparison with the footprint vs RNA-sequencing samples, we need to again ensure that these labels for this "Type" factor are listed in the correct alphabetical order to ensure we are performing a footprint *VS* RNA-sequencing comparison to reflect translation efficiency.

.. ident with TABs
.. code-block:: python

  **sample_info.txt**

    Sample    Condition   Type
    s1_fp     a_WT        RPF
    s1_rna    a_WT        RNA
    s2_fp     a_WT        RPF
    s2_rna    a_WT        RNA
    s3_fp     b_EXP       RPF
    s3_rna    b_EXP       RNA
    s4_fp     b_EXP       RPF
    s4_rna    b_EXP       RNA

| The alphabetical order of the factor names (i.e., "Condition", "Type") does not matter. Instead, according to the DESeq2 documentation, these design factors are `evaluated in the order listed <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset>`_. However, changes to the order will cause negligible differences in output. For example, if we scramble the order the factors are listed in the design formula, we obtain essentially the same output:

.. ident with TABs
.. code-block:: python

  $ xpresspipe diffxpress -i tm1.tsv -s tm_deseq.txt \
                          --design Type+Condition+Type:Condition

          "baseMean"           "log2FoldChange"    "lfcSE"                "stat"              "pvalue"                "padj"
  "ATF4"  3283.07267363348     2.5427843106451     0.134284452518271      18.9358057687216    5.78241665951195e-80    5.02954601044349e-76


  $ xpresspipe diffxpress -i tm_counts.tsv -s tm_deseq.txt \
                          --design Condition+Type+Condition:Type

          "baseMean"           "log2FoldChange"    "lfcSE"                "stat"              "pvalue"                "padj"
  "ATF4"  3283.07267363348     2.54278431064905    0.134284452463494      18.9358057764753    5.78241580816713e-80    5.02954526994377e-76



| For more information on factor levels and design parameters, please see the `DESeq2 documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs>`_ and `this note <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#can-i-run-deseq2-to-contrast-the-levels-of-many-groups>`_. Any standard design formula that will work in DESeq2 will work within the XPRESSpipe wrapper, as long as the formatted described above is followed.

| Other possible variations to DESeq2 analysis are available `here <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variations-to-the-standard-workflow>`_, but not all will be compatible with the XPRESSpipe wrapper. In general, the XPRESSpipe wrapper is best suited to simple multi-factor design (Experimental vs Wild-type, Footprints vs RNA-sequencing, plus any other factors relevant to your experiment). For advice in preparing your design formula differently than in the examples listed below, please reach out to us `here <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe diffxpress --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path/filename.tsv\>`, :data:`--input \<path/filename.tsv\>`
     - Path and file name of expression counts matrix
   * - :data:`-s \<path/filename.tsv\>`, :data:`--sample \<path/filename.tsv\>`
     - Path and file name of sample information matrix
   * - :data:`--design \<formula\>`
     - Design formula for differential expression analysis (spaces in command line are conserved in input string. DO NOT INCLUDE ~ OR SPACES IN FORMULA IN COMMAND LINE, will be automatically added)

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`--suppress_version_check`
    - Suppress version checks and other features that require internet access during processing
  * - :data:`--shrink`
    - Provide argument to perform shrinkage of effect size on log fold changes. Useful for visualization and ranking of hits


--------------------------------------------
Example 1 -- Analyze ribosome profiling data
--------------------------------------------
| The source files can be found `here <https://github.com/XPRESSyourself/xpressyourself_manuscript/tree/main/isrib_analysis/isrib_de/xpresspipe_data_deseq2>`_.
| If we want to perform differential expression of translation efficiency for ribosome profiling data, we need to provide :data:`Condition` and :data:`Type` factor columns in the :data:`sample_info` file. If we want to include the :data:`RPF` / :data:`RNA` comparison to account for translation efficiency, we would need to include these factor label as a column to ensure the appropriate :data:`RPF` / :data:`RNA` evaluation.

.. ident with TABs
.. code-block:: python

  **tm_counts.tsv**

          ribo_untr_a  ribo_untr_b  ribo_tm_a  ribo_tm_b  untr_a_hek  untr_b_hek  tm_a_hek  tm_b_hek
  A1BG    29           43           21         11         67          73          56        85
  A2M     3            5            2          2          73          57          32        37
  AAAS    1441         1981         934        601        1144        1067        1012      1124
  AACS    575          727          310        192        351         335         220       291
  AADAT   98           120          51         29         322         315         192       292


  **tm_deseq.txt**

  Sample         Condition       Type
  untr_a_hek     UNTR            RNA
  untr_b_hek     UNTR            RNA
  ribo_untr_a    UNTR            RPF
  ribo_untr_b    UNTR            RPF
  tm_a_hek       TM              RNA
  tm_b_hek       TM              RNA
  ribo_tm_a      TM              RPF
  ribo_tm_b      TM              RPF

.. code-block:: shell

  $ xpresspipe diffxpress -i counts_data.tsv --sample sample_info.txt --design Type+Condition+Type:Condition

| The output of this analysis will perform differential expression that reflects both :data:`TM` vs :data:`UNTR` *and* :data:`RPF` (footprints) vs :data:`RNA`.


.. ident with TABs
.. code-block:: python

  **tm_counts_diffx.tsv**

          baseMean	 log2FoldChange	       lfcSE             stat	         pvalue	         padj
  ATF4	  3283.072674	 2.542784311	       0.134284453	 18.93580577	 5.78E-80	 5.03E-76
  PTP4A1  460.6444433	 2.473962772	       0.185061193	 13.36834986	 9.26E-41	 4.03E-37
  SPEN	  7902.554413	 1.192124338	       0.109445545	 10.89239713	 1.25E-27	 3.63E-24
  RPS15A  1823.967865	 -1.391099082	       0.152069954	 -9.147757652	 5.81E-20	 1.26E-16
  DYNC1H1 11985.60418	 0.85282198	       0.094425503	  9.031691164	 1.69E-19	 2.56E-16

| From this output, we can focus on the :data:`log2FoldChange` and :data:`padj` columns. From this output, we see that ATF4 is the most significantly upregulated gene by translation efficiency between the TM and UNTR conditions, which is what we expect (see the `XPRESSyourself manuscript <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007625>`_ for further discussion of this example). Further explanations of the other columns of this output can be found in the `DESeq2 documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.


---------------------------------
Example 2 -- Analyze RNA-seq data
---------------------------------
| For a standard two-condition RNA-seq experiment comparison, we are only interested in the differential expression of :data:`EXP` vs :data:`WT`. To ensure this comparison if performed correctly, we need to force these :data:`Condition` factor labels to be alphabetical. We will thus rename them :data:`b_EXP` and :data:`a_WT` and do the following:

.. ident with TABs
.. code-block:: python

  **expression_counts.tsv**

                  s1  s2  s3  s4  ...
  ENSG00000227232 66  59  1   82  ...
  ENSG00000240361 35  0   7   72  ...
  ENSG00000238009 20  70  85  78  ...
  ENSG00000241860 96  7   93  38  ...
  ENSG00000187634 73  41  92  77  ...


  **sample_info.tsv**

  Sample  Condition
  s1      a_WT
  s2      a_WT
  s3      a_WT
  s4      a_WT
  s5      b_EXP
  s6      b_EXP
  s7      b_EXP
  s8      b_EXP

.. code-block:: shell

  $ xpresspipe diffxpress -i test_r/test_dataset.tsv --sample test_r/sample_info.tsv --design Condition


-------------------------------------------------------------------------
Example 3 -- Analyze RNA-seq data that was prepared in different batches
-------------------------------------------------------------------------
| If samples were performed in multiple batches and you would like to control for batch effect, you can add a :data:`Batch` factor column and provide different batch labels. This example below will control for batch effect and compare :data:`EXP` vs :data:`WT` expression.
| See the `DESeq2 documentation example <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start>`_ for further information.

.. ident with TABs
.. code-block:: python

  **expression_counts.tsv**

                  s1  s2  s3  s4  ...
  ENSG00000227232 66  59  1   82  ...
  ENSG00000240361 35  0   7   72  ...
  ENSG00000238009 20  70  85  78  ...
  ENSG00000241860 96  7   93  38  ...
  ENSG00000187634 73  41  92  77  ...


  **sample_info.tsv**

  Sample  Condition Batch
  s1      a_WT      batch1
  s2      a_WT      batch1
  s3      a_WT      batch1
  s4      a_WT      batch1
  s5      b_EXP     batch2
  s6      b_EXP     batch2
  s7      b_EXP     batch2
  s8      b_EXP     batch2

.. code-block:: shell

  $ xpresspipe diffxpress -i test_r/test_dataset.tsv --sample test_r/sample_info.tsv --design Batch+Condition




======================
rRNA Probe
======================
| Ribosome RNA (rRNA) contamination is common in RNA-seq library preparation. As the bulk of RNA in a cell at any given time is dedicated to rRNA, and as these rRNA sequences are relatively few and therefore highly repeated, depletion of these sequences is often desired in order to have better depth of coverage of non-rRNA sequences. In order to facilitate this depletion, many commercial kits are available that target specific rRNA sequences for depletion, or that enrich mRNA polyA tails. However, and especially in the case of ribosome profiling experiments, where RNA is digested to create ribosome footprints that commercial depletion kits won't detect and polyA selection kits are inoperable as footprints will not have the requisite polyA sequence. To this end, `custom rRNA probes <https://www.ncbi.nlm.nih.gov/pubmed/28579404>`_ are recommended, and the :data:`rrnaProbe` sub-module was designed to facilitate this process.
| :data:`rrnaProbe` works by doing the following:
| 1. Run FASTQC to detect over-represented sequences
| 2. Collate these sequences to determine consensus fragments
| 3. Output rank ordered list of over-represented fragments within the appropriate length range to target for depletion
| NOTE: BLAST capability to verify over-represented consensus fragments are indeed rRNA sequences is not yet incorporated, so any sequences that will be used as probes should be BLAST-verified first.

.. code-block:: shell

  $ xpresspipe rrnaProbe --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to zipped FASTQC files
   * - :data:`-o \</path/filename\>, --output \</path/filename\>`
     - Path and file name to write output

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--suppress_version_check`
     - Suppress version checks and other features that require internet access during processing
   * - :data:`-m \<value\>, --min_overlap \<value\>`
     - Minimum number of bases that must match on a side to combine sequences (default: 5)
   * - :data:`--footprint_only`
     - Only take zip files that are ribosome profiling footprints (file names must contain "FP", "RPF", or "FOOTPRINT")

-----------
Examples
-----------
| **Example 1 -- Generate rank-ordered list of over-represented sequences**

.. ident with TABs
.. code-block:: python

  $ xpresspipe rrnaProbe -i riboprof_out/fastqc_out/ -o riboprof_out/sequences.txt --footprint_only

  TTGATGATTCATAATAACTTTTCGAATCGCAT    514832
  TATAAATCATTTGTATACGACTTAGAT         121739
  TTGATGATTCATAATAACTTTTCGAATCGCAT    15776
  TTTGATGATTCATAATAACTTTTCGAATCGCAC   33325
  ATAAATCATTTGTATACGACTTAGAC          13603
