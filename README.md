# ![XPRESSpipe](https://raw.githubusercontent.com/XPRESSyourself/XPRESSpipe/master/docs/content/xpresspipe.png)


### A toolkit for navigating and analyzing gene expression datasets

[![Build Status](https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master)](https://travis-ci.org/XPRESSyourself/XPRESSpipe)
[![Documentation Status](https://readthedocs.org/projects/xpresspipe/badge/?version=latest)](https://xpresspipe.readthedocs.io/en/latest/?badge=latest)
[![PyPi Status](https://img.shields.io/pypi/v/XPRESSpipe.svg)](https://pypi.org/project/XPRESSpipe/)

##### Find documentation [here](https://xpresspipe.readthedocs.io/en/latest/)

-----

### Development Notes:
<b><i>XPRESSpipe is still in alpha production</i></b>   

### Citation:    
```
Berg, JA (2019). XPRESSyourself suite: Gene expression processing and analysis made easy. https://github.com/XPRESSyourself.
```

### Installation:   
Installation options not currently available   
```
pip install xpresspipe
```
```
conda install -c bioconda xpresspipe
```

### QuickStart:   
```
$ xpresspipe riboprof -i /path/to/raw/data/ -o /path/to/output/ -r /path/to/reference/ ...
```

### Important Notes:    
In order for many of the XPRESSpipe functions to perform properly and for the output to be reliable after alignment (except for generation of a raw counts table), recommended file naming conventions must be followed.

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
