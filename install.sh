#!/bin/bash

git clone https://github.com/XPRESSyourself/XPRESSpipe.git

conda install -f XPRESSpipe/requirements.yml

cd XPRESSpipe; python setup.py install; cd ../

pip uninstall xpresstools

Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("Rsubread")'
Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("dupRadar")'
Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("DESeq2")'

git clone https://github.com/XPRESSyourself/XPRESStools.git; cd XPRESStools; python setup.py install
