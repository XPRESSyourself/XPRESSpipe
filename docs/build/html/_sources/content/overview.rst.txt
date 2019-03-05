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
Batch Effect Normalization
===========================
| **Batch Effect**: When using larger RNAseq libraries where multiple people prepared sequence libraries or where libraries were prepared on different days, batch effect normalization must be performed to remove these effects on sequence counts and make libraries that are comparable to one another. This step can be performed using the :data:`xpresspipe batchNormalize` sub-module.
| **Sample Variability**: Other sample inherent biases, such as variable total read counts per sample on the same sequencing chip, should be normalized using the :data:`xpresspipe normalizeMatrix` sub-module.
