FROM continuumio/miniconda3
RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

RUN git clone https://github.com/XPRESSyourself/XPRESSpipe.git

RUN conda env create -f ./XPRESSpipe/requirements.yml
RUN echo "source activate $(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)/bin:$PATH

RUN cd XPRESSpipe; python setup.py install

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev r-base \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("Rsubread")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("dupRadar")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("DESeq2")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("GenomicFeatures")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("GenomicAlignments")'
RUN Rscript -e 'install.packages("data.table", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'library(devtools)' -e 'install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)'

RUN git clone https://github.com/XPRESSyourself/XPRESStools.git; cd XPRESStools; python setup.py install

ENTRYPOINT ["xpresspipe"]
