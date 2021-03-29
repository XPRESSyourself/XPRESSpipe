.. _install_link:

############
Installation
############

=================================
Install XPRESSpipe
=================================
| 1. Open your command line interface and install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_, if not already installed.

.. code-block:: shell

  # If on a MacOS
  $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

  # If on a LinuxOS
  $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  $ bash ~/Miniconda3-latest-MacOSX-x86_64.sh

  # Enter yes for all successive prompts and allow the script to install Conda into your path
  # After installation, the install script can be removed
  rm ~/Miniconda3-latest-MacOSX-x86_64.sh


| 2. Download the latest version of XPRESSpipe by executing the lines of code in the code block below. Replace the URL for the version of XPRESSpipe for whatever version you want (these can be found under the :data:`releases` tab on the `XPRESSpipe GitHub repository <https://github.com/XPRESSyourself/XPRESSpipe/releases>`_).

.. code-block:: shell

  $ cd ~
  $ curl -L -O https://github.com/XPRESSyourself/XPRESSpipe/archive/v0.6.3.tar.gz
  $ tar xvzf v0.6.3.tar.gz
  $ cd XPRESSpipe-0.6.3

| 3. Install XPRESSpipe dependencies via Conda and activate the XPRESSpipe environment:

.. code-block:: shell

  $ conda env create --name xpresspipe -f requirements.yml
  $ conda activate xpresspipe

| 4. This installation method will create a separate environment for XPRESSpipe and all its dependencies to live in. Each time you open the command line, you will need to type :data:`conda activate xpresspipe` to use XPRESSpipe
| 5. Install XPRESSpipe and test that the installation was successful:

.. note::
  :data:`v0.6.3` and later employs the :data:`bash install.sh` method for installing XPRESSpipe. If using :data:`v0.6.2` or earlier, you should instead run :data:`pip install .`

.. code-block:: shell

  $ bash install.sh
  $ xpresspipe test

| If a summary menu appeared in the command line interface, it means we are good to go! Congrats! You are almost ready to use XPRESSpipe!
|
| You can run :data:`xpresspipe --help` to see a list of the available modules within XPRESSpipe. To see specific parameters for a module, type :data:`xpresspipe <module_name> --help`.


==============================================================
Install in a supercomputing environment
==============================================================
| If the dependencies for XPRESSpipe were installed to a conda environment as above, you will need to add a couple lines to you bash script to submit the supercomputing job.
| For example, if using a SLURM job scheduler, you should include the following after the :data:`#SBATCH` lines and before any calls to XPRESSpipe in the slurm script, as below:

.. code-block:: shell

  #!/bin/bash
  #SBATCH --time=72:00:00
  #SBATCH --nodes=1
  #SBATCH ...

  source $(conda info --base)/etc/profile.d/conda.sh
  source activate xpresspipe

  ... rest of the script
