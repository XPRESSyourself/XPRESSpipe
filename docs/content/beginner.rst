############
Beginners
############

=================================
First Steps
=================================
| If this is your first time doing any programming, congratulations! You are embarking upon a very rewarding path. As with learning any new natural language, there is a learning curve associated with learning a computer language. While XPRESSpipe is aimed at reducing majority of the overhead associated with processing this data, using this software will still require some effort, just as would learning any new language or laboratory technique.

| XPRESSpipe is used through something called the `command line interface <https://en.wikipedia.org/wiki/Command-line_interface>`_ (or CLI), or what some people refer to as `"The Matrix" <https://www.youtube.com/watch?v=kqUR3KtWbTk>`_. This may seem daunting, but luckily, several free online courses are available to quickly catch you up to speed on some of the basics that will be required to use this software. We recommend Codecademy's CLI course, which you can find `here <https://www.codecademy.com/learn/learn-the-command-line>`_ and should take only a couple of hours (Codecademy estimates ~10 hours, but you probably don't need to finish the course to use XPRESSpipe. The purpose of this is to help you become more comfortable with the command line).

| Once, you're ready to jump into the command line, we can get rolling! For the steps below, we're going to assume we are on an Mac operating system and provide examples under this pretext, but this software is compatible with any Linux-like operating system and the syntax is largely the same (sorry Windows users!).

=================================
Install XPRESSpipe
=================================
| - Let's enter the command line.
| 1. Click on the Finder icon the top right side of the screen on your Mac (or wherever else it might be located)
| 2. Type "Terminal" into the search bar and click on the app icon

| - Great! Now we are in the command line interface. As a review, anything followed by a "$" in the command line is a command and you can execute each command by pressing Enter after typing. You can also auto-complete file names using Tab. But be careful, **space and characters must be typed exactly and commands are case-sensitive**.
| - First, let's make sure the required dependency is installed. If the :data:`pip` command does not work, then read `this <https://pip.pypa.io/en/stable/installing/>`_ for more information.

.. code-block:: shell

  $ pip install setuptools

| - Let's get the latest version of XPRESSpipe by executing the lines of code in the code block below. Replace the URL for the version of XPRESSpipe for whatever version you want (these can be found under the :data:`releases` tab on the XPRESSpipe GitHub repository).

.. code-block:: shell

  $ cd ~
  $ curl -O https://github.com/XPRESSyourself/XPRESSpipe/archive/XPRESSpipe-v0.1.3b2.tar.gz
  $ tar xvzf XPRESSpipe-v0.1.3b2.tar.gz
  $ cd XPRESSpipe-v0.1.3b2
  $ python setup.py install

| - Let's test that this worked by executing the following:

.. code-block:: shell

  $ xpresspipe -h

| - If a help menu appeared in the command line interface, it means we are good to go! Congrats! You are almost ready to use XPRESSpipe!


=================================
Generate Reference Files
=================================
| - Before we can actually use XPRESSpipe to process our raw RNA-seq data, we need to create a reference directory. Directory is just programmer lingo for a folder. In this example, we will be working with human-derived RNA-seq data, so let's perform the following in the command line:

.. code-block:: shell

  $ cd ~/Desktop
  $ mkdir reference_folder
  $ mkdir reference_folder/fasta_files

| 1. The first command helped us navigate to the Desktop. The "~" icon is a shortcut for the User directory, and every directory needs to be separated by a "/"
| 2. The second command created a new folder in the Desktop directory called :data:`reference_folder`
| 3. The third command created a new folder in the reference directory for intermediate reference files

| - Now let's get the reference files. We're going to do this directly in the command line, but if you have trouble with this, I will explain an alternative afterwards. Quick note, because the next lines of code are a bit long, I used the "\" character to indicate I am continuing the command in the next line. You do not need this in executing the command, they just help make the code a little more readable.

.. code-block:: shell

  $ cd reference_folder/
  $ curl ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz -o transcripts.gtf.gz
  $ gzip -d *.gz
  $ cd fasta_files/
  $ for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; \
      do curl -O ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz; \
      done
  $ gzip -d *.gz
  $ cd ../

| 1. We navigated into the reference folder, downloaded a GTF reference file and unzipped it, then navigated to the :data:`fasta_file` directory to download the raw reference data and unzipped it. Finally, we returned to the main reference directory.
| 2. If this didn't work, we can navigate to `Ensembl <https://www.ensembl.org/>`_ to get the relevant data. We need to get the `GTF file <ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz>`_ and `each chromosome sequence file <ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/>`_. You can follow the links to download these files and then move them into your reference folder. The link to the chromosome sequence files actually contains more files than we need. We just need the files that start with :data:`Homo_sapiens.GRCh38.dna.chromosome`. If these files were zipped with a :data:`.zip` or :data:`.gz` extension, double click each file to unzip them.

| - Now we need to curate these references files into something the sequencing alignment software can use. Since we are using ribosome profiling data, we want a reference that will allow us to `avoid mapping to the 5' and 3' ends of genes <https://www.cell.com/cms/10.1016/j.celrep.2016.01.043/attachment/257faf34-ff8f-4071-a642-bfdb531c75b8/mmc1>`_. We also don't want to align to anything but protein coding genes. Finally, we want to quantify to the longest transcript. This last bit just helps the software avoid confusion when a gene has multiple splice variants to choose from. Since this is short read sequencing, we also want to factor this into the curation of the reference (see the :data:`--sjdbOverhang` argument below).

.. code-block:: shell

  $ docker run jordanberg/xpresspipe curateReference --output ./ \
                                                      --fasta fasta_files/ \
                                                      --gtf ./transcripts.gtf \
                                                      --longest_transcript \
                                                      --protein_coding \
                                                      --truncate \
                                                      --sjdbOverhang 49

| - The truncation option is only necessary when using XPRESSpipe to process ribosome profiling samples and their associated RNA-seq samples.
| - If interested in quantifying miRNA, etc, leave out the :data:`--protein_coding` argument.
| - If running sequencing where the read (single-end) or mates not equal to 100 bp, you will want to change the :data:`--sjdbOverhang` argument to be the length of one of the paired-end reads - 1, so if we ran 2x100bp sequencing, we would specify :data:`--sjdbOverhang 99`
| - This may take awhile, and as we will discuss later, you may want to run these steps on a supercomputer, but this will serve as a preliminary guide for now.

=================================
Process Raw Sequencing Files
=================================
| - Now let's get our raw data. Let's follow the following instructions:
| 1. Make a new folder, something called :data:`raw_data` or whatever you like and place your data there.
| 2. Make sure the files follow proper naming conventions (see naming conventions :ref:`here <general-usage>`)
| 3. Now let's process the data
| 4. Let's also create a folder called something like :data:`output`
| 5. Also, make sure you have the 3' adaptor sequence handy used when generating your sequencing library
| 6. We'll feed the program the new GTF file that contains only longest transcript, protein coding, truncated references generating in the reference curation step
| 7. We'll give the experiment a name and also specify what `method of sample normalization <https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>`_ we want performed on the count data
| 8. We also need to specify the :data:`--sjdbOverhang` amount we fed into the reference curation step, so in this case we will use :data:`--sjdbOverhang 49`

.. code-block:: shell

  $ docker run jordanberg/xpresspipe riboprof --input raw_data/ \
                                              --output output/ \
                                              --reference reference_folder/ \
                                              --gtf reference_folder/transcripts_longestTranscript_proteinCoding_truncated.gtf
                                              --experiment riboprof_test
                                              --adaptor CTGTAGGCACCATCAAT
                                              --method RPKM
                                              --sjdbOverhang 49

| - If you are running a lot of files, especially for human samples, this may take a lot of time. We recommend running this on some kind of server. A situation like yeast with few samples may be feasible to run on a personal computer, but will likely also take some time

======================
Explore the Data
======================
| - Once the data is finished processing, we can start exploring the output.

------------------
Sequencing Metrics
------------------


------------------
Library Complexity
------------------


-------------------
Metagene Analysis
-------------------


--------------------------------
Periodicity (Ribosome Profiling)
--------------------------------


----------------------------------
Count Data and Downstream Analysis
----------------------------------












=======================
Supercomputing
=======================

---------------
Getting Started
---------------


---------------
Load XPRESSpipe
---------------





---------------
Load Data
---------------


----------------
Curate Reference
----------------




---------------
Process Data
---------------



--------------
Retrieve Data
--------------
