FROM continuumio/miniconda3
RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

RUN git clone https://github.com/XPRESSyourself/XPRESStools.git; cd XPRESStools; python setup.py install

RUN git clone https://github.com/XPRESSyourself/XPRESSpipe.git

RUN conda env create -f ./XPRESSpipe/requirements.yml
RUN echo "source activate $(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /XPRESSpipe/requirements.yml | cut -d' ' -f2)/bin:$PATH

RUN cd XPRESSpipe; python setup.py install

ENTRYPOINT ["xpresspipe"]
