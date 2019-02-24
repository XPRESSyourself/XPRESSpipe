############
Installation
############

===================
MacOS Installation
===================

See local_install.sh in the `resources <https://github.com/RiboPipe/ribopipe/tree/master/resources>`_ folder for interactive script.

1)  RiboPipe requires use of command line. Execute the following lines of code in `Terminal <https://www.imore.com/how-use-terminal-mac-when-you-have-no-idea-where-start>`_ (on Mac, open Spotlight and type 'Terminal'):

2)  Install python3, wget, git, and git-lfs if not already done. You can check in command line by typing in the name of the package and checking if your system recognizes the package name.

.. code-block:: shell

  $ sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
  $ echo "export PATH='$(brew --prefix)/bin:$(brew --prefix)/sbin'":'"$PATH"' >> ~/.bash_profile
  $ brew install python3 wget git git-lfs
  $ git lfs install

3)  You may need to manually point your system to python3. You can check this by typing :data:`python` in the command line and seeing if it is running python3. At this time, RiboPipe only works when Python3 is set as the default. If it is not, do the following:

.. code-block:: shell

  $ echo "alias python="python3" >> ~/.bash_profile

4)  Download `Conda <https://www.anaconda.com/download/#macos>`_, a package manager, for your operating system. Double click the `.pkg` file if on MacOS, the `.exe` file on Windows, or follow these `instructions <https://conda.io/docs/user-guide/install/linux.html#install-linux-silent>`_ on Linux.

5)  Download RiboPipe dependencies via conda:

.. code-block:: shell

  $ conda install -y -c bioconda setuptools fastqc fastx_toolkit htseq picard samtools hisat2 star bedtools deeptools scipy plastid pandas numpy matplotlib seaborn pysam=0.14

6)  Install RiboPipe dependencies via pip (should have come pre-installed with python3.4 or greater):

.. code-block:: shell

  $ pip install multiqc

7)  Download RiboPipe:

.. code-block:: shell

  $ git clone https://github.com/RiboPipe/ribopipe.git
  $ cd ribopipe
  $ cd python setup.py install

8)  Or, to download specific version:

.. code-block:: shell

  $ tag='v0.1.4-beta'
  $ wget https://github.com/RiboPipe/ribopipe/archive/$tag.zip
  $ unzip ribopipe-${tag:1}.zip
  $ mv ribopipe-${tag:1} ribopipe
  $ cd ribopipe
  $ cd python setup.py install

9)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing ribopipe script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/RiboPipe-0.1.5b0-py3.6.egg
  Processing dependencies for RiboPipe==0.1.5b0
  Finished processing dependencies for RiboPipe==0.1.5b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile

10) Test installation:

.. code-block:: shell

  $ ribopipe --help


=================
HPC Installation
=================

See hpc_install.sh in the `resources <https://github.com/RiboPipe/ribopipe/tree/master/resources>`_ folder for interactive script. While the resource manager can install these for you, we will show you how to manually install all dependencies. These instructions may vary slightly from HPC to HPC.

1)  Remove all pre-loaded software:

.. code-block:: shell

  $ module purge

2)  Install brew and related dependencies:
Install python3, wget, git, and git-lfs if not already done. You can check in command line by typing in the name of the package and checking if your system recognizes the package name.

.. code-block:: shell

  $ sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
  $ echo "export PATH='$(brew --prefix)/bin:$(brew --prefix)/sbin'":'"$PATH"' >> ~/.bash_profile
  $ brew install python3 wget git git-lfs
  $ git lfs install

3)  You may need to manually point your system to python3. You can check this by typing :data:`python` in the command line and seeing if it is running python3. At this time, RiboPipe only works when Python3 is set as the default. If it is not, do the following:

.. code-block:: shell

  $ echo "alias python="python3" >> ~/.bash_profile

4)  Install anaconda (instructions retrieves most recent version as of time of writing):

.. code-block:: shell

  $ wget https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh
  $ chmod 700 Anaconda3-5.3.0-Linux-x86_64.sh
  $ ./Anaconda3-5.3.0-Linux-x86_64.sh -b -p $HOME/.local/bin -s
  $ export PATH="/uufs/chpc.utah.edu/common/home/$USER/.local/bin:$PATH"
  $ conda install -y -c bioconda setuptools fastqc fastx_toolkit htseq picard samtools star bedtools deeptools scipy plastid pandas numpy matplotlib seaborn pysam=0.14

5)  Install pip related dependencies (should have come pre-installed with python3.4 or greater):

.. code-block:: shell

  $ pip install multiqc

6)  To download current repository:

.. code-block:: shell

  $ git clone https://github.com/RiboPipe/ribopipe.git
  $ cd ribopipe
  $ python setup.py install --prefix ~/.local

7)  Or, to download specific version

.. code-block:: shell

  $ tag='v0.1.4-beta'
  $ wget https://github.com/RiboPipe/ribopipe/archive/$tag.zip
  $ unzip ribopipe-${tag:1}.zip
  $ mv ribopipe-${tag:1} ribopipe
  $ cd ribopipe
  $ cd python setup.py install --prefix ~/.local

8)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing ribopipe script to /uufs/chpc.utah.edu/common/home/$USER/.local/bin

  Installed /uufs/chpc.utah.edu/common/home/$USER/.local/lib/python3.5/site-packages/RiboPipe-0.1.5b0-py3.5.egg
  Processing dependencies for RiboPipe==0.1.5b0
  Finished processing dependencies for RiboPipe==0.1.5b0

  $ echo "export PATH='/uufs/chpc.utah.edu/common/home/$USER/.local/bin:$PATH' >> ~/.bash_profile

9) Test installation:

.. code-block:: shell

  $ ribopipe --help


==================
Test Installation
==================

To test installation, run the following command:

.. code-block:: shell

  $ ribopipe riboseq -i /path/to/ribopipe/test/ -o /path/you/create/ -r yeast -e ingolia_2015 \
  -p STAR -a CTGTAGGCACCATCAAT --platform ILLUMINA --count_cutoff 32 \
  -s a_wild-type_DED1_replicate_1_15_deg c_ded1-cs_replicate_1_15_deg

If no errors are produced by the output, the installation was successful.


================
Additional Help
================

Manually installing a package:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sometimes a package may need to be manually installed. In these cases, a pattern as follows may be used (example given is loading on HPC).

.. code-block:: shell

  $ wget https://github.com/simon-anders/htseq/archive/release_0.11.0.zip
  $ unzip htseq-release_0.11.0.zip
  $ rm htseq-release_0.11.0.zip
  $ cd htseq-release_0.11.0
  $ python setup.py install --prefix ~/.local
  $ cd ../
  $ echo "export PATH='/uufs/chpc.utah.edu/common/home/$USER/.local/bin:$PATH' >> ~/.bash_profile"

Or...

.. code-block:: shell

  $ git clone https://github.com/simon-anders/htseq.git
  $ cd htseq
  $ python setup.py install --prefix ~/.local
  $ cd ../
  $ echo "export PATH='/uufs/chpc.utah.edu/common/home/$USER/.local/bin:$PATH' >> ~/.bash_profile"

Getting publicly available raw data from GEO:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Raw data from previous studies that have been made publicly available can be accessed through the `GEO database <https://www.ncbi.nlm.nih.gov/geo/>`_

Please see `this example script <https://github.com/RiboPipe/ribopipe/blob/master/resources/run_riboseq_GSE66411_test.sh>`_ for examples of how to retrieve this data. This `thread <https://www.biostars.org/p/111040/>`_ is also helpful.
