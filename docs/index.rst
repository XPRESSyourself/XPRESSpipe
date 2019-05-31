##############
XPRESSpipe
##############
|build-status| |docs| |Docker|

=================
About
=================
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is a part of the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite of sequencing tools. `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is an automated, efficient, and flexible pipeline for processing RNA-seq raw data, performing quality control analysis on the data, and preparing data for further downstream analysis.  `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is currently capable of handling single-end (SE), paired-end (PE), and ribosome profiling data. Features include the ability to trim, align, and count sequence reads.
| Other useful features include:
| - Curate reference files for alignment
| - Prepare transcriptome reference file to contain only protein coding transcripts, longest isoform only, and/or truncate transcripts
| - Identify over-abundant sequences in the read data to design possible rRNA depletion probes
| - Perform transcript meta-analysis, ribosome profiling periodicity analysis, and library complexity analysis
| - Output aligned reads in a variety of formats for downstream applications
|
| Other analyses can be performed by `XPRESSplot <https://github.com/XPRESSyourself/XPRESSplot>`_. Please read the relevant documentation for more information.
|
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ and the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite is developed and maintained by Jordan Berg in the `Rutter Lab <https://biochem.utah.edu/rutter/index.html>`_ @ the `University of Utah <https://www.utah.edu/>`_, along with other collaborators.
|
| If you have limited or no computational experience, please see our :ref:`beginners_link`

=================
Table of contents
=================
.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 1

   content/overview
   content/beginner
   content/installation
   content/general-usage
   content/quickstart
   content/reference-building
   content/trim
   content/align
   content/count
   content/formatting
   content/normalization
   content/analysis
   content/seRNAseq
   content/peRNAseq
   content/riboseq
   content/faqs

=======
License
=======
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is perpetually open access under a GNU General Public License (v3.0).

==========
Questions?
==========
| If you have questions, requests, or bugs to report, please use the `XPRESSpipe issues forum <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.



.. |build-status| image:: https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master
    :target: https://travis-ci.org/XPRESSyourself/XPRESSpipe
    :alt: Build Status

.. |codecov| image:: https://codecov.io/gh/XPRESSyourself/XPRESSpipe/XPRESSpipe.svg?branch=master
    :target: https://codecov.io/gh/XPRESSyourself/XPRESSpipe
    :alt: Code Coverage

.. |docs| image:: https://readthedocs.org/projects/xpresspipe/badge/?version=latest
    :target: https://xpresspipe.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Docker| image:: https://img.shields.io/static/v1.svg?label=docker&message=dowload&color=informational
    :target: https://cloud.docker.com/repository/docker/jordanberg/xpresspipe/general
    :alt: Docker
