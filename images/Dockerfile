FROM continuumio/miniconda3

RUN apt-get update \
 && apt-get -y install --no-install-recommends \
   build-essential \
   libxml2-dev \
   zlib1g-dev \
   gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN echo "Starting clone of repo..."
RUN git clone https://github.com/XPRESSyourself/XPRESSpipe.git
RUN conda init bash
RUN conda update conda
RUN conda env create -n xpresspipe python=3.7 -f ./XPRESSpipe/requirements.yml
ENV PATH /opt/conda/envs/xpresspipe/bin:$PATH
RUN echo "source activate xpresspipe" > ~/.bashrc
ENV CONDA_DEFAULT_ENV xpresspipe
ENV PATH /opt/conda/envs/$(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)/bin:$PATH
RUN cd XPRESSpipe; python setup.py install

RUN apt-get update \
 && apt-get -y install --no-install-recommends \
   build-essential \
   libxml2-dev \
   zlib1g-dev \
   gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("DESeq2", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("Rsubread", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("dupRadar", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("GenomicFeatures", dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install("GenomicAlignments", dependencies=TRUE)'
RUN Rscript -e 'install.packages("data.table", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("BiocParallel")'
RUN Rscript -e 'BiocManager::install("Rhtslib")'
RUN Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'library(devtools)' -e 'devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)'

ENTRYPOINT ["xpresspipe"]
