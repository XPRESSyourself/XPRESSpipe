############
Overview
############

====================
SE and PE RNA-seq
====================
| The XPRESSpipe pipeline is flexibly designed to be able to process and perform preliminary analyses on single-end (SE) or paired-end RNA-seq sequence read. Raw data is most often generated in the form of a `.fastq <http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm>`_ or .txt file. This data is useful in determining the transcriptional landscape of a population of cells or of single cells. Other qualities, such as microRNA abundance, splice events, and sequence variants can also be detected and analyzed.

====================
Ribosome Profiling
====================
| `Ribosome profiling <https://en.wikipedia.org/wiki/Ribosome_profiling>`_ utilizes Next Generation Sequencing (NGS) to provide a detailed picture of the protein translation landscape within cells. Cells are lysed, translating ribosomes are isolated, and the ribosome protected mRNA fragments (ribosome footprints are integrated into a SE RNA-seq library. The library is then sequenced and processed similarly to a single-end RNA-seq run, with some exceptions:

| - **5' and 3' ribosome footprint bias**: 5' and 3' footprint bias has been well-documented in any ribosome profiling library preparation method, arising from the longer ribosome initiation and termination steps. Therefore, it is advised that when quantifying reads to a reference transcriptome, one avoids mapping to the first 45-50 nucleotides of exon 1 of each transcript. By using the :data:`xpresspipe truncate` sub-module, one can create a truncated transcriptome reference file (GTF) for use when quantifying their aligned reads.

| - **rRNA contamination**: As pervasive ribosomal RNA (rRNA) can contaminate footprint libraries due to most commercial kits' inability to deplete these various fragmented rRNAs inherent in the ribosome profiling methodology, it is advised to create depletion probes for dominant rRNA fragment species in ribosome profiling libraries for a given organism. By using the :data:`xpresspipe rrnaProbe` sub-module, one can collate these dominant rRNA species and create depletion probes to prevent their incorporation into future sequence libraries.

| See this `paper <https://www.ncbi.nlm.nih.gov/pubmed/28579404>`_ for a recent discussion and detailed protocol of the technique.

.. image:: riboseq_overview.png
   :width: 600
   :align: center

===========================
Software
===========================
| XPRESSpipe aims to use the most up-to-date iterations of a software type required for high-throughput genomics processing and analysis. In order to design XPRESSpipe, we referred to a variety of benchmarking studies comparing the software we chose with others to verify that performance and speed would be optimal. Below is a rationale for each package chosen. As software continues to improve and benchmarking studies are published, XPRESSpipe will also be updated to reflect these improvements.

| `fastp <https://github.com/OpenGene/fastp>`_ -- Read pre-processing
| While external benchmarking has not been performed to our knowledge in recent years on read pre-processing tools, we chose to use fastp as it is fast, and (at least from self-reports) has reliable output. We also favored fastp as it is able to handle more recent trends in RNA-seq, such as trimming of `unique molecular identifiers (UMIs) <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4933-1>`_.

| `STAR <https://github.com/alexdobin/STAR>`_ -- Masking and Alignment
| A recent `benchmarking paper <https://www.nature.com/articles/nmeth.4106>`_ showed that STAR outperformed other comparable tools in speed and performance, increasing the number of correctly aligned reads, while reducing the number of falsely called reads as is the case with several other packages. Additionally, in the XPRESS pipeline, we have included an option to perform a preliminary masking step, where one can remove sequences that often confound mRNA expression levels (i.e. rRNAs, tRNAs).

| `samtools <https://github.com/samtools/samtools>`_/`bedtools <https://github.com/arq5x/bedtools2>`_/`deepTools <https://github.com/deeptools/deepTools>`_ -- Alignment file handling
| These tools handle the alignment file processing before quantification to identify PCR duplicates (optional), remove non-uniquely aligned reads, and so on.

| `Cufflinks <https://github.com/cole-trapnell-lab/cufflinks>`_ -- Read quantification (Recommended)
| A recent `benchmarking paper <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0734-x>`_ showed evidence that Cufflinks using default parameters performed the best compared to other read quantification tools.
| XPRESSpipe uses Cufflinks v2.1.1 as Cufflinks v2.2.1 suffers from a persistent Seg Fault 11 error on MacOS. No significant changes effecting quantification have occurred between these versions.

| `HTSeq <https://github.com/simon-anders/htseq>`_ -- Read quantification
| While not performing optimally compared to other tools, HTSeq has been included for archival purposes as the current `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ pipeline uses HTSeq. This is also a included in cases where the user requires raw counts as quantification output.

| `dupRadar <https://bioconductor.org/packages/release/bioc/html/dupRadar.html>`_ -- Library Complexity
| dupRadar is a stable, easy to use tool for estimating library size complexity and doesn't suffer from systematic software issues like other tools that contain a similar module

| `SVA <http://bioconductor.org/packages/release/bioc/html/sva.html>`_ -- Known Library Batch Correction
| Used for correcting for known batch effects between samples (i.e. samples prepared on different days, by different people, etc.)

| `DESeq2 <http://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ -- Differential Expression Analysis
| Perform differential expression analysis on the data.

| `MultiQC <https://github.com/ewels/MultiQC>`_ -- Summary reports
| MultiQC gathers log output from fastp, STAR, and Cufflinks to provide the user with a easy to view summary of their processed data at each step


=======================
Methodology
=======================

| Additionally, we seek to provide the best methodology for high-throughput sequencing processing, and explain key components below.

| **Transcriptomic Reference Files**
| Read quantification often requires a transcriptome reference file in order to know what alignment coordinates map to what genes. We introduce a suite of GTF modification tools included in XPRESSpipe that we will briefly discuss:
| - Isoforms: GTF files contain records for every isoform of a gene. However, since these isoforms all contain overlapping regions, many quantification tools count a read mapping to one of these regions as a multi-mapper and either penalizes it or discards it completely. A common way to handle this is by taking only the longest transcript for each gene during quantification. This can be performed with :data:`xpresspipe modifyGTF -l`.
| - Protein Coding: When calculating mRNA expression levels, sample normalization to reduce technical bias from RNA-seq platforms is important. However, highly-abundant rRNAs can confound these metrics. Therefore, we provide an option to create a GTF file with only protein-coding annotated genes as input for quantification using :data:`xpresspipe modifyGTF -p`.
| - Ribosome Profiling Bias: During translation, there are three steps: 1) Initiation, 2) Elongation, and 3) Termination. There is usually a pause during Initiation and Termination, which will present itself as systematic spikes on the 5' and 3' ends of each transcript for ribosome profiling reads. A way to correct for the kinetics of initiation and termination and measure translational capacity itself is to `avoid mapping reads to the first 15 codons and last 5 codons of a transcript <https://www.ncbi.nlm.nih.gov/pubmed/28579404>`_. :data:`xpresspipe modifyGTF -t` handles this by searching the exon space of each transcript and pruning the given amounts off of each so that these regions are considered non-coding space. This process is performed recursively, so that if you were trimming 45 nt from the 5' end and exon 1 was only 30 nt, exon 1 would be removed and exon 2 would be trimmed by 15 nt.

| **PCR De-Duplication**
| During sequence library creation, a PCR amplification step is common in order to produce enough sequence material, but often, different reads are amplified differentially.  When UMIs are not used, these duplication events can artificially model higher expression of a transcript that had favorable amplification conditions. We therefore include an optional PCR de-duplication step for experiments not using UMIs. Be warned, this can introduce `additional biases <https://www.ncbi.nlm.nih.gov/pubmed/30001700>`_ and should be used with caution. Performing library complexity analysis on the samples should indicate whether or not computational de-duplication should be performed.

| **Meta-Analysis**
| - Periodicity: A helpful metric of ribosome profiling libraries is looking at the characteristic 3 nt/1 codon stepping of the translating ribosome.

| - Metagene: Metagene analysis takes the read coverage across all transcripts in a sample and compiles their distribution along a representative transcript. This is useful in identifying any systematic 5' or 3' biases in the library preparation step.

| - Read distribution: Once reads are trimmed of low quality bases or adaptor sequences, looking at the distribution of read lengths can be helpful in identifying that the expected RNA was incorporated into the library. This is especially useful in ribosome profiling datasets, where ideally all reads isolated and incorporated into the library should fall within the 21-30 nt range.
