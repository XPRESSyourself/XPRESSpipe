FROM continuumio/miniconda3
RUN conda create -n env python=3.7
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
RUN conda install -c r r-devtools
RUN git clone https://github.com/XPRESSyourself/XPRESSpipe.git

RUN conda env create -f ./XPRESSpipe/requirements.yml
RUN echo "source activate $(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)/bin:$PATH

RUN cd XPRESSpipe; python setup.py install

RUN apt-get update \
 && apt-get -y install --no-install-recommends \
   libfftw3-dev \
   r-base \
   build-essential \
   libcurl4-gnutls-dev \
   libxml2-dev \
   gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos = "http://cran.us.r-project.org")}'
RUN Rscript -e 'BiocManager::install("DESeq2", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("Rsubread", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("dupRadar", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("GenomicFeatures", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("GenomicAlignments", dependencies=TRUE)'
RUN Rscript -e 'install.packages("data.table", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("BiocParallel")'
RUN Rscript -e 'library(devtools)' -e 'devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)'


RUN git clone https://github.com/XPRESSyourself/XPRESStools.git; cd XPRESStools; python setup.py install

ENTRYPOINT ["xpresspipe"]
