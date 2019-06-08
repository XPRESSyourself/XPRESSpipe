############
Installation
############

=================================
Install from Source (Recommended)
=================================
| The majority of installation is automated within the setup script.
| 1. Assuming you have a more recent version of Python, you should have :data:`pip` installed. If not, please follow these `instructions <https://pip.pypa.io/en/stable/installing/>`_.
| 2. Install the required :data:`setuptools` package that will automate setup.

.. code-block:: shell

  $ pip install setuptools

| 3. Download XPRESSpipe and execute the XPRESSpipe installation script.

.. code-block:: shell

  $ cd ~
  $ curl -O https://github.com/XPRESSyourself/XPRESSpipe/archive/XPRESSpipe-v0.1.3b2.tar.gz
  $ tar xvzf XPRESSpipe-v0.1.3b2.tar.gz
  $ cd XPRESSpipe-v0.1.3b2
  $ python setup.py install

| This will install XPRESSpipe as well as all required dependencies. If `Anaconda <https://www.anaconda.com/distribution/>`_ has not been installed by the user, this will take care of installation. Anaconda is a package manager and will install the other software packages XPRESSpipe relies on.

| 4. Test installation by executing the following in the command line:

.. code-block:: shell

  $ xpresspipe -h

| 5. If your system can't find the :data:`xpresspipe` command, try the following. At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing xpresspipe script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/xpresspipe-0.0.1b0-py3.6.egg
  Processing dependencies for XPRESSpipe==0.0.1b0
  Finished processing dependencies for XPRESSpipe==0.0.1b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile
  $ echo "export PYTHONPATH='/Users/$USERNAME/anaconda3/lib/python3.6/site-packages'"

======================
Docker Container
======================
| NOTE: The Docker containerization is still under development.
| XPRESSpipe is available as a fully independent `Docker <https://www.docker.com/>`_ container. By downloading the Docker software and using the below command, a image of XPRESSpipe and all associated dependencies localized to a single file for ease of use.

| 1. Download XPRESSpipe Docker image:

.. code-block:: shell

  $ docker image pull jordanberg/xpresspipe:latest

| 2. Run XPRESSpipe:

.. code-block:: shell

  $ docker run jordanberg/xpresspipe --help

| If the help menu prints, XPRESSpipe if functioning properly and you can replace the :data:`--help` option with the appropriate sub-module and arguments.



====================
Dependencies
====================
In both installation options specified above, the installation of these packages should be handled automatically. If not, please make sure you have `Anaconda <https://www.anaconda.com/distribution/>`_ installed and try the installation again if using the Installation from Source option. If there is an issue with the Docker distribution, please submit an `issue <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.
| - fastp
| - STAR
| - samtools
| - bedtools
| - deepTools
| - FastQC
| - Cufflinks
| - HTSeq
| - pandas
| - numpy
| - biopython
| - scipy
| - plotly
| - R
| - dupRadar
| - DESeq2
| - matplotlib
| - XPRESSplot
| - MultiQC
