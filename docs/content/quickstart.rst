##########
Quickstart
##########

To do a full install on a local machine or HPC, see the `Installation <installation.html>`_ page.


Singularity
^^^^^^^^^^^
For a faster way to get started, consider using a Singularity image for RiboPipe, hosted `here <https://github.com/RiboPipe/ribopipe-singularity>`_.

1)  `Singularity <https://www.sylabs.io/docs/>`_ offers a fast, reproducible way of running RiboPipe where all dependencies and OS are bundled into a single disk image. These singularity "containers" are widely used in cloud computing.
`Download <https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps>`_ if you do not already have it (must be on a Linux OS). Often, cloud computing environments come with this software pre-installed.

2)  Download a RiboPipe singularity container:

.. code-block:: shell

  $ singularity pull library://sylabsed/linux/ribopipe

3)  Run Ribopipe:

.. code-block:: shell

  $ raw_data=/path/to/raw/data
  $ output_data=/path/to/output/data
  $ singularity exec ribopipe.sif riboseq -i $raw_data -o $output_data ...


Local
^^^^^
1) Move raw data to directory of choice
2) Create empty output directory
3) Run ribopipe:

.. code-block:: shell

  $ raw_data=/path/to/raw/data
  $ output_data=/path/to/output/data
  $ ribopipe riboseq -i $raw_data -o $output_data ...

4) Collect raw_counts.csv output in $output_data/assembly/counts and edit `sample_info.csv <https://github.com/RiboPipe/ribopipe/blob/master/resources/diffex_template.csv>`_
5) Run diffex:

.. code-block:: shell

  $ ribopipe_path=/path/to/ribopipe
  $ ribopipe diffex -i $output_data/assembly/counts/raw_counts.csv -d $ribopipe_path/resources/sample_info.csv -o output_name --type riboseq


HPC
^^^
1) Modify hpc_run_template.sh in the `resources <https://github.com/RiboPipe/ribopipe/tree/master/resources>`_
folder for an example script for submitting the pipeline job to the HPC and make sure dependencies
listed in this script are on the HPC system, else they need to be locally installed
2)  Run the script by executing the following:

.. code-block:: shell

  $ sbatch hpc_run_template.sh

If you want the slurm output file to be sent to the SLURM directory to avoid storage space issues on your interactive node, then in the `#SBATCH -o slurmjob-%j` line, replace it with the path to your SLURM directory:

.. code-block:: shell

  $ #SBATCH -o /scratch/general/lustre/INPUT_USER_ID_HERE/slurmjob-%j
