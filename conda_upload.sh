#!/bin/bash

#Inspired by Qiusheng Wu, https://wetlands.io

PKG_NAME='xpresspipe'
USER='j-berg'
# adjust the Python versions you would like to build
array=( 2.7 3.5 3.6 3.7 )
echo "Building conda package ..."
cd ~
conda skeleton pypi $PKG_NAME
cd $PKG_NAME
cd ~
# building conda packages
for i in "${array[@]}"
do
	conda-build --python $i $PKG_NAME
done
# convert package to other platforms
cd ~
platforms=( osx-64 linux-32 linux-64 win-32 win-64 )
find $HOME/conda-bld/linux-64/ -name *.tar.bz2 | while read file
do
    echo $file
    #conda convert --platform all $file  -o $HOME/conda-bld/
    for platform in "${platforms[@]}"
    do
       conda convert --platform $platform $file  -o $HOME/conda-bld/
    done
done
# upload packages to conda
find $HOME/conda-bld/ -name *.tar.bz2 | while read file
do
    echo $file
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l nightly $file --force
done
echo "Building conda package done!"
