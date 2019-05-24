# ![XPRESSpipe](https://raw.githubusercontent.com/XPRESSyourself/XPRESSpipe/master/docs/content/xpresspipe.png)


### An alignment and analysis pipeline for RNAseq data

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESSpipe)
[![codecov.io](https://codecov.io/gh/XPRESSyourself/XPRESSpipe/XPRESSpipe.svg?branch=master)](https://codecov.io/gh/XPRESSyourself/XPRESSpipe)
[![Documentation Status](https://readthedocs.org/projects/xpresspipe/badge/?version=latest)](https://xpresspipe.readthedocs.io/en/latest/?badge=latest)
[![Docker](https://img.shields.io/static/v1.svg?label=docker&message=dowload&color=informational)](https://cloud.docker.com/repository/docker/jordanberg/xpresspipe/general)

-----
Please refer to the [documentation](https://xpresspipe.readthedocs.io/en/latest/?badge=latest) for more in depth details. 

### Development Notes:
- <b><i>XPRESSpipe is still in beta production</i></b>  
- The current release [XPRESSpipe-v0.1.3b2](https://github.com/XPRESSyourself/XPRESSpipe/releases/tag/XPRESSpipe-v0.1.3b2) is running relatively stable on MacOS and Linux (including HPCs)
  - The meta-analysis plotting seems to work on a local MacOS, but not when running on a Linux HPC, will hopefully have kinks worked out in the next couple of weeks
  - Yet to incorporate UMI handling, representative gene housekeeping, along with some other features in the near future

### Citation:    
```
Berg, JA, et. al. (2019). XPRESSyourself: Automating and Democratizing High-Throughput Sequencing. https://github.com/XPRESSyourself.
```

### Installation:   
#### Installing from source
1. Installation requires Python (distributed with most operating systems automatically) and setuptools. If you have a more current version of Python, you can install setuptools as follows:
```
$ pip install setuptools
```
If this does not work, please refer to this [site](https://pip.pypa.io/en/stable/installing/) for more information
2. Get XPRESSpipe by downloading and unpacking the most recent archive found [here](https://github.com/XPRESSyourself/XPRESSpipe/releases)
3. Unzip the folder and navigate to the appropriate directory in the command line
```
$ cd /path/to/XPRESSpipe
```
4. Install XPRESSpipe
```
$ python setup.py install
```
5. Test the installation
```
$ xpresspipe -h
```
If the help menu is not displayed, try adding the path where you installed XPRESSpipe to the system PATH
```
$ echo 'export PATH=$PATH:/path/to/xpresspipe' >> ~/.bash_profile
```
If you do not have a file names `~/.bash_profile`, try looking for one called `~/.profile`


#### Using a Docker container
1. [Install Docker](https://docs.docker.com/v17.12/install/)
2. Download the XPRESSpipe Docker container
```
$ docker pull docker push jordanberg/xpresspipe:latest
```
3. Run the Docker container
```
docker run jordanberg/xpresspipe --help
```

### QuickStart:   
```
$ xpresspipe riboprof -i /path/to/raw/data/ -o /path/to/output/ -r /path/to/reference/ ...
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
