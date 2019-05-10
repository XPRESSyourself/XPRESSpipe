#!/bin/bash

git clone https://github.com/XPRESSyourself/XPRESSpipe.git

conda install --file XPRESSpipe/requirements.yml

cd XPRESSpipe; python setup.py install; cd ../

pip uninstall -y xpresstools

git clone https://github.com/XPRESSyourself/XPRESStools.git; cd XPRESStools; python setup.py install


Rscript install_dependencies.r
