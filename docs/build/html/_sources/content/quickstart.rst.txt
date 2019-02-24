##########
Quickstart
##########

1) `Install <installation>`_ RiboPipe.

Local
^^^^^

2) Move raw data to directory of choice
3) Create output directory
4) Run ribopipe:

.. code-block:: shell

  $ raw_data=/path/to/raw/data
  $ output_data=/path/to/output/data
  $ ribopipe riboseq -i $raw_data -o $output_data

5) Collect raw_counts.csv output in $output_data/assembly/counts and edit `sample_info.csv <https://github.com/j-berg/ribopipe/blob/master/resources/diffex_template.csv>`_
6) Run diffex:

.. code-block:: shell

  $ ribopipe_path=/path/to/ribopipe
  $ ribopipe diffex -i $output_data/assembly/counts/raw_counts.csv -d $ribopipe_path/resources/sample_info.csv -o output_name --type riboseq


HPC
^^^

2) Modify hpc_run_template.sh in the `resources <https://github.com/j-berg/ribopipe/resources/>`_
folder for an example script for submitting the pipeline job to the HPC and make sure dependencies
listed in this script are on the HPC system, else they need to be locally installed
3)  Run the script by executing the following:

.. code-block:: shell

  $ sbatch hpc_run_template.sh

If you want the slurm output file to be sent to the SLURM directory to avoid storage space issues on your interactive node, then in the `#SBATCH -o slurmjob-%j` line, replace it with the path to your SLURM directory:

.. code-block:: shell

  $ #SBATCH -o /scratch/general/lustre/INPUT_USER_ID_HERE/slurmjob-%j
