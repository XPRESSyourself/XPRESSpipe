#!/bin/bash
# Adapted from mikecormier (https://github.com/gogetdata/ggd-cli/blob/master/.circleci/setup.sh)
# The MIT License (MIT)

# Copyright (c) 2016 gogetdata

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set -exo pipefail

WORKSPACE=$(pwd)

# Set path
echo "export PATH=$WORKSPACE/anaconda/bin:$PATH" >> $BASH_ENV
source $BASH_ENV

## Passed from .circleci/config.yml (Only 2 or 3 permited)
pythonversion=$1
if (( $pythonversion != 2 && $pythonversion != 3 ))
then
    echo -e "\nERROR: Python 2 or 3 designation required. Python version $pythonversion was supplied. Please correct and run again\n"
    exit 1
fi

# setup conda and dependencies
if [[ ! -d $WORKSPACE/anaconda ]]; then
    mkdir -p $WORKSPACE


    # step 1: download and install anaconda
    if [[ $OSTYPE == darwin* ]]; then
        tag="MacOSX"
        tag2="darwin"
    elif [[ $OSTYPE == linux* ]]; then
        tag="Linux"
        tag2="linux"
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi

    curl -O https://repo.continuum.io/miniconda/Miniconda$pythonversion-latest-$tag-x86_64.sh
    sudo bash Miniconda$pythonversion-latest-$tag-x86_64.sh -b -p $WORKSPACE/anaconda/
    sudo chown -R $USER $WORKSPACE/anaconda/
    mkdir -p $WORKSPACE/anaconda/conda-bld/$tag-64

    # step 2: setup channels
    conda config --system --add channels defaults
    conda config --system --add channels conda-forge
    conda config --system --add channels bioconda
    conda config --system --add channels r

    # step 3: install xpresspipe requirements
    conda config --set always_yes yes --set changeps1 no --set show_channel_urls true
    conda info -a
    conda env create -n test --file requirements.yml
    source activate test
    conda info --envs
    conda list
    pip install coverage

    # step 5: cleanup
    conda clean -y --all

    # Add local channel as highest priority
    mkdir -p $WORKSPACE/miniconda/conda-bld/{noarch,linux-64,osx-64}
    conda index $WORKSPACE/miniconda/conda-bld
    conda config --system --add channels file://$WORKSPACE/miniconda/conda-bld

fi

conda config --get

ls $WORKSPACE/miniconda/conda-bld
ls $WORKSPACE/miniconda/conda-bld/noarch
