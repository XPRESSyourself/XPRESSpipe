#!/bin/bash

OS=$(uname -s)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "Installing XPRESSpipe for ${OS}..."

# Install XPRESSpipe
pip install $DIR

# Install fastp-lite
cd fastp_lite
if [ ${OS} == "Darwin" ]
then
    make -f Makefile_macOS
else
    make -f Makefile_Linux
fi
cd fastp_lite ../xpresspipe
cd $DIR

# Install cufflinks
if [ ${OS} == "Darwin" ]
then
  curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.OSX_x86_64.tar.gz -o $DIR/cufflinks.tar.gz
  tar -zxvf $DIR/cufflinks.tar.gz
  rm $DIR/cufflinks.tar.gz;
  mv $DIR/cufflinks-2.1.1.OSX_x86_64 $DIR/cufflinks
  mv $DIR/cufflinks/cufflinks $DIR/xpresspipe
  echo "Cufflinks installed";
else
  curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz -o $DIR/cufflinks.tar.gz
  tar -zxvf $DIR/cufflinks.tar.gz
  rm $DIR/cufflinks.tar.gz;
  mv $DIR/cufflinks-2.1.1.Linux_x86_64 $DIR/cufflinks
  mv $DIR/cufflinks/cufflinks $DIR/xpresspipe
  echo "Cufflinks installed";
fi

echo "XPRESSpipe installation complete."
