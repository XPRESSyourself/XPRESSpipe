# ![XPRESSpipe](https://raw.githubusercontent.com/XPRESSyourself/XPRESSpipe/master/docs/content/xpresspipe.png)


### An alignment and analysis pipeline for RNAseq data

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESSpipe)
[![Documentation Status](https://readthedocs.org/projects/xpresspipe/badge/?version=latest)](https://xpresspipe.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/170939943.svg)](https://zenodo.org/badge/latestdoi/170939943)
<!--[![Anaconda-Server Badge](https://anaconda.org/bioconda/xpresspipe/badges/version.svg)](https://anaconda.org/bioconda/xpresspipe)-->
<!--[![codecov.io](https://codecov.io/gh/XPRESSyourself/XPRESSpipe/XPRESSpipe.svg?branch=master)](https://codecov.io/gh/XPRESSyourself/XPRESSpipe)-->
<!--[![Docker](https://img.shields.io/static/v1.svg?label=docker&message=dowload&color=informational)](https://cloud.docker.com/repository/docker/jordanberg/xpresspipe/general)-->

-----
Please refer to the [documentation](https://xpresspipe.readthedocs.io/en/latest/?badge=latest) for more in depth details.

### Citation:    
```
Berg JA, et. al. (2020). XPRESSyourself: Enhancing, standardizing, and
automating ribosome profiling computational analyses yields improved insight
into data. PLoS Comp Biol. doi: https://doi.org/10.1371/journal.pcbi.1007625
```

### Installation:   
#### Installing from source
The following is a short tutorial showing you how to install XPRESSpipe:   
[![asciicast](https://asciinema.org/a/262192.svg)](https://asciinema.org/a/262192?speed=4)

- Make sure you let Anaconda set up the PATH info for you.
- If the help menu is not displayed when testing, try adding the path where you installed XPRESSpipe to the system PATH
```
$ echo 'export PATH=$PATH:/path/to/xpresspipe' >> ~/.bash_profile
```
- If you do not have a file names `~/.bash_profile`, try looking for one called `~/.profile`
- The commands used in the video above are summarized here:
```
$ curl -L -O https://github.com/XPRESSyourself/XPRESSpipe/archive/v0.2.3b0.zip
$ unzip v0.2.3b0.zip
$ cd XPRESSpipe-0.2.3b0/
$ conda env create -f requirements.yml
$ conda activate xpresspipe
$ python setup.py install
$ xpresspipe -h
$ xpresspipe test
```
- Be sure to specify the correct release version in the first URL

### QuickStart:   
- Reference building   
[![asciicast](https://asciinema.org/a/256340.svg)](https://asciinema.org/a/256340?speed=4)

- Running XPRESSpipe on sequence data   
[![asciicast](https://asciinema.org/a/256343.svg)](https://asciinema.org/a/256343?speed=4)

- You can also use the XPRESSpipe command builder and executor for reference curation or running the pipeline by executing the following:
```
$ xpresspipe build
```


### Important Notes:    
#### Basic Starting Input
- `input` directory with raw sequence data
  - Sequence data files should be `FASTQ` format and end in `.fastq` or `.fq` and can be `.zip` or `.gz` compressed
- An empty `output` directory
- A `reference` directory (see documentation for `curateReference` for more details)

#### Naming Conventions
In order for ordered output after alignment (except for generation of a raw counts table), recommended file naming conventions should be followed.

1. Download your raw sequence data and place in a folder -- this folder should contain all the sequence data and nothing else.
2. Make sure files follow a pattern naming scheme. For example, if you had 3 genetic backgrounds of ribosome profiling data, the naming scheme would go as follows:
```
ExperimentName_BackgroundA_FP.fastq(.qz)
ExperimentName_BackgroundA_RNA.fastq(.qz)
ExperimentName_BackgroundB_FP.fastq(.qz)
ExperimentName_BackgroundB_RNA.fastq(.qz)
ExperimentName_BackgroundC_FP.fastq(.qz)
ExperimentName_BackgroundC_RNA.fastq(.qz)
```
3. If the sample names are replicates, their sample number needs to be indicated.
4. If you want the final count table to be in a particular order and the samples ordered that way are not alphabetically, append a letter in front of the sample name to force this ordering.
```
ExperimentName_a_WT.fastq(.qz)
ExperimentName_a_WT.fastq(.qz)
ExperimentName_b_exType.fastq(.qz)
ExperimentName_b_exType.fastq(.qz)
```
5. If you have replicates:
```
ExperimentName_a_WT_1.fastq(.qz)
ExperimentName_a_WT_1.fastq(.qz)
ExperimentName_a_WT_2.fastq(.qz)
ExperimentName_a_WT_2.fastq(.qz)
ExperimentName_b_exType_1.fastq(.qz)
ExperimentName_b_exType_1.fastq(.qz)
ExperimentName_b_exType_2.fastq(.qz)
ExperimentName_b_exType_2.fastq(.qz)
```


### Running a test dataset:
- We can run a test dataset as in the [associated manuscript](https://www.biorxiv.org/content/10.1101/704320v1) by downloading the FASTQ files from [GSE65778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65778) using the [SRAtoolkit](https://www.biostars.org/p/111040/).
- We can curate the reference like so:
```
$ xpresspipe curateReference -o /path/to/reference -f /path/to/reference/genome_fastas -g /path/to/reference/transcripts.gtf -p -t --sjdbOverhang 49
```
- And we can process the dataset like so:
```
xpresspipe riboseq -i /path/to/input -o /path/to/output -r /path/to/reference/ --gtf /path/to/reference//transcripts_CT.gtf -e isrib_test_study -a CTGTAGGCACCATCAAT --sjdbOverhang 49
```
- The above steps will be very computationally intensive, so we recommend running this on a supercomputing cluster
- Scripts used to analyze this data can be found [here](https://github.com/XPRESSyourself/xpressyourself_manuscript/blob/master/isrib_analysis/isrib_analysis.py) and [here](https://github.com/XPRESSyourself/xpressyourself_manuscript/blob/master/isrib_analysis/isrib_de/run_de.sh) and [here](https://github.com/XPRESSyourself/xpressyourself_manuscript/blob/master/isrib_analysis/isrib_de/isrib_de_analysis.py)

- Alternatively, smaller test datasets can be found within the XPRESSpipe `tests` folder and an outline of commands to run can be found [here](https://github.com/XPRESSyourself/XPRESSpipe/blob/master/tests/test_pipelines.py)


### Updates
Information on updates to the software can be found [here](https://github.com/XPRESSyourself/XPRESSpipe/releases).
