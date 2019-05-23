############
Overview
############

====================
SE and PE RNAseq
====================
| The XPRESSpipe pipeline is flexibly designed to be able to process and perform preliminary analyses on single-end (SE) or paired-end RNAseq sequence read. Raw data is most often generated in the form of a `.fastq <http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm>`_ or .txt file. This data is useful in determining the transcriptional landscape of a population of cells or of single cells. Other qualities, such as microRNA abundance, splice events, and sequence variants can also be detected and analyzed.

====================
Ribosome Profiling
====================
| `Ribosome profiling <https://en.wikipedia.org/wiki/Ribosome_profiling>`_ utilizes Next Generation Sequencing (NGS) to provide a detailed picture of the protein translation landscape within cells. Cells are lysed, translating ribosomes are isolated, and the ribosome protected mRNA fragments (ribosome footprints are integrated into a SE RNAseq library. The library is then sequenced and processed similarly to a single-end RNAseq run, with some exceptions:

| - **5' ribosome footprint bias**: 5' footprint bias has been well-documented in any sequence library preparation method. Therefore, it is advised that when counting reads to a reference transcriptome, one avoids mapping to the first 45-50 nucleotides of exon 1 of each transcript. By using the :data:`xpresspipe truncate` sub-module, one can create truncated transcriptome reference files (GTF) for counting their aligned reads.

| - **rRNA contamination**: As pervasive ribosomal RNA (rRNA) can contaminate footprint libraries due to most commercial kits' inability to deplete these various fragmented rRNAs inherent in the ribosome profiling methodology, it is advised to create depletion probes that are dominant in ribosome profiling libraries for a given organism. By using the :data:`xpresspipe rrnaProbe` sub-module, one can collate these dominant rRNA species and create depletion probes for these fragment species.

| See this `paper <https://www.ncbi.nlm.nih.gov/pubmed/28579404>`_ for a recent discussion and detailed protocol of the technique.

.. image:: riboseq_overview.png
   :width: 600
   :align: center

===========================
Software
===========================
| XPRESSpipe aims to use the most up-to-date iterations of a software type required for high-throughput genomics processing and analysis. In order to design XPRESSpipe, we referred to a variety of benchmarking studies comparing the software we chose with others to verify that performance and speed would be optimal. Below is a rationale for each package chosen. As software continues to improve and benchmarking studies are published, XPRESSpipe will also be updated to reflect these improvements.

| `fastp <https://github.com/OpenGene/fastp>`_ -- Read pre-processing

| `STAR <https://github.com/alexdobin/STAR>`_ -- Masking and Alignment

| `samtools <https://github.com/samtools/samtools>`_/`bedtools <https://github.com/arq5x/bedtools2>`_/`deepTools <https://github.com/deeptools/deepTools>`_ -- Alignment file handling

| `Cufflinks <https://github.com/cole-trapnell-lab/cufflinks>`_ -- Read quantification (Recommended)
| Using v2.1.1 as v2.2.1 suffers from a persistent Seg Fault 11 error on MacOS. No significant changes occurred between these versions influencing performance of Cufflinks quantification

| `HTSeq <https://github.com/simon-anders/htseq>`_ -- Read quantification

| `dupRadar <https://bioconductor.org/packages/release/bioc/html/dupRadar.html>`_ -- Library complexity

| `MultiQC <https://github.com/ewels/MultiQC>`_ -- Summary reports






=======================
Methodology
=======================

| Additionally, we seek to provide the best methodology for high-throughput sequencing processing, and explain key components below.

| GTF mods

| Dedup


| Metagene


| Periodicity

| Read distribution
