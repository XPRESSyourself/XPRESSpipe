############
Installation
############

=================================
Install XPRESSpipe
=================================
| - Let's enter the command line.
| 1. Click on the Finder icon the top right side of the screen on your Mac (or wherever else it might be located)
| 2. Type "Terminal" into the search bar and click on the app icon

| - Great! Now we are in the command line interface. As a review, anything followed by a "$" in the command line is a command and you can execute each command by pressing Enter after typing. You can also auto-complete file names using Tab. But be careful, **space and characters must be typed exactly and commands are case-sensitive**.
| - First, let's install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_, a package manager that will help us download all required dependencies.

.. code-block:: shell

  # If on a MacOS
  $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

  # If on a LinuxOS
  $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  $ bash ~/Miniconda3-latest-MacOSX-x86_64.sh

  # Enter yes for all successive prompts and allow the script to install Conda into your path
  # After installation, the install script can be removed
  rm ~/Miniconda3-latest-MacOSX-x86_64.sh



| - Let's get the latest version of XPRESSpipe by executing the lines of code in the code block below. Replace the URL for the version of XPRESSpipe for whatever version you want (these can be found under the :data:`releases` tab on the `XPRESSpipe GitHub repository <https://github.com/XPRESSyourself/XPRESSpipe/releases>`_).

.. code-block:: shell

  $ cd ~
  $ curl -L -O https://github.com/XPRESSyourself/XPRESSpipe/archive/v0.2.1b0.tar.gz
  $ tar xvzf v0.2.1b0.tar.gz
  $ cd XPRESSpipe-0.2.1b0

| - Now we can give Conda the dependency file to process for us from XPRESSpipe:

.. code-block:: shell

  $ conda env create --name xpresspipe -f requirements.yml
  # conda activate xpresspipe


| - This installation method will create a separate environment for XPRESSpipe and all its dependencies to live in. Each time you open the command line, you will need to type :data:`conda activate xpresspipe` to use XPRESSpipe
| - Let's install XPRESSpipe and test that the installation was successful by executing the following:

.. code-block:: shell

  $ python setup.py install
  $ xpresspipe test

| - If a summary menu appeared in the command line interface, it means we are good to go! Congrats! You are almost ready to use XPRESSpipe!
| - You can run :data:`xpresspipe --help` to see a list of the available modules within XPRESSpipe. To see specific parameters for a module, type :data:`xpresspipe <module_name> --help`.


------------------------------
Important!
------------------------------
| Currently, R is struggling to auto-install riboWaltz for p-site QC. If you run into issues in the time-being, please open an R session in your environment and run the following commands:

.. code-block:: R

  install.packages("devtools", repos = "http://cran.us.r-project.org")

  devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = c("Depends", "Imports", "LinkingTo"))

| And then run the following to test installation:

.. code-block:: R

  library(riboWaltz)


==============================================================
Using XPRESSpipe in a supercomputing environment (i.e. SLURM)
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
