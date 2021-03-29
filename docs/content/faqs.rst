############################
FAQs
############################

| If you have questions, requests, or bugs to report, please use the `XPRESSpipe issues forum <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.

| As this software is used more, we will compile some of the most commonly asked questions here...

-----------------------------------------------------------
A step of the pipeline is erroring for no apparent reason
-----------------------------------------------------------
| First, please check the output in your terminal, along with in the log file. If the step that the pipeline breaks on does not output any useful information, check that the required dependencies were installed correctly. For example, when we were testing the the :data:`geneCoverage` module on a supercomputing cluster, the pipeline responded saying it couldn't find the appropriate index file. It turned out the R package, GenomicFeatures was not downloaded due to issues with the rtracklayers package. For this situation, we fixed it by uninstalling Anaconda and reinstalling the dependencies, as below:

.. code-block:: shell

    # Run each of these steps. If a command doesn't work, skip to the next one
    $ conda install anaconda-clean
    $ anaconda-clean --yes
    $ rm -rf ~/miniconda
    $ rm ~/.condarc
    $ rm -r ~/.conda/

.. code-block:: shell

    $ cd ~
    $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash ~/Miniconda3-latest-Linux-x86_64.sh
    $ conda env base --file XPRESSpipe/requirements.yml

| During Anaconda installation, reply :data:`yes` to all prompts. If you wish to install the XPRESSpipe dependencies to their own environment, replace :data:`base` with :data:`your_environment_name_here` in the last step. If XPRESSpipe continues to malfunction after completion of these steps, please reach out to us on the `XPRESSpipe issues forum <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.

---------------------------------------------------------------------------
The pipeline breaks because of a segmentation fault during alignment.
---------------------------------------------------------------------------
| Occasionally, depending on allocation of CPUs, 32 virtual CPUs may be available, but only 16 are configured. This may lead to memory overloads by trying to use more than configured, as the large index files will be temporarily copied to each processing core. If this is the case, provide the :data:`max_processors` with the number in the log file stated as available. For a computing node with 64 GB of RAM available, we generally see that 20 CPUs is stable. See log example below:

.. code-block:: shell

    sh: line 1: 70311 Segmentation fault      STAR --runThreadN 30 ...

    or

    WARNING: fastp uses up to 16 threads although you specified 32
