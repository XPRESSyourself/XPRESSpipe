############
Installation
############

====================
PyPi Install
====================
| This feature is not yet available...
|
| 1)  Install XPRESSpipe and associated dependencies via pip:

.. code-block:: shell

  $ pip install xpresspipe

====================
Conda Install
====================
| This feature is not yet available...
|
| 1)  Install XPRESSpipe and associated dependencies via conda:

.. code-block:: shell

  $ conda install -y -c bioconda xpresspipe

======================
Manual install
======================

| 1)  Or download XPRESSpipe manually:

.. code-block:: shell

  $ git clone https://github.com/XPRESSyourself/xpresspipe.git
  $ cd XPRESSpipe
  $ cd python setup.py install

| 2)  Or, to download specific version:

.. code-block:: shell

  $ tag='v0.0.1-beta'
  $ wget https://github.com/XPRESSyourself/XPRESSpipe/archive/$tag.zip
  $ unzip xpresspipe-${tag:1}.zip
  $ mv xpresspipe-${tag:1} xpresspipe
  $ cd xpresspipe
  $ cd python setup.py install

| 3)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing xpresspipe script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/xpresspipe-0.0.1b0-py3.6.egg
  Processing dependencies for XPRESSpipe==0.0.1b0
  Finished processing dependencies for XPRESSpipe==0.0.1b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile

====================
Dependencies
====================
| - FastP
| - STAR
| - samtools
| - bedtools
| - deeptools
| - xpresspipe
| - picard
| - fastqc
| - plastid
| - gtfToGenePred
| - htseq
