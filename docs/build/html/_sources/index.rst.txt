##############
XPRESSpipe
##############
|build-status| |docs| |pypi| |doi|

=================
About
=================
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is a part of the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite of sequencing tools. `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is an automated and flexible pipeline for processing RNAseq raw data and preparing for sequence analysis. `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is capable of handling single-end (SE), paired-end (PE), and ribosome profiling SE data. Features include the ability to trim, align, and count sequence reads. Other useful features include:
| - Create a 5' truncated transcript reference file to counting reads (useful in ribosome profiling to avoid quantifying the persistent 5' ribosome footprint bias)
| - Automated creation of `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ alignment reference files
| - Identification of over-abundant sequences in library prep for depletion probe design of those sequences and library quality improvement
| - Meta-analyses, such as metagene profiles for each processed file and 3-codon periodicity of most-abundant footprint length (useful in ribosome profiling library quality control)
|
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ and the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite is developed and maintained by Jordan Berg in the `Rutter Lab <https://biochem.utah.edu/rutter/index.html>`_ @ the `University of Utah <https://www.utah.edu/>`_, along with other collaborators.

=================
Table of contents
=================
.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 1

   content/overview
   content/installation
   content/general-usage
   content/quickstart
   content/reference-building
   content/seRNAseq
   content/peRNAseq
   content/riboprof
   content/trim
   content/align
   content/count
   content/quality
   content/formatting
   content/normalization

=======
License
=======
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is freely available under a GNU General Public License (v3.0).

==========
Questions?
==========
| If you have questions, requests, or bugs to report, please use the `Github issues forum <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.



.. |build-status| image:: https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master
    :target: https://travis-ci.org/XPRESSyourself/XPRESSpipe
    :alt: Build Status

.. |docs| image:: https://readthedocs.org/projects/xpresspipe/badge/?version=latest
    :target: https://xpresspipe.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |pypi| image:: https://img.shields.io/pypi/v/XPRESSpipe.svg
    :target: https://pypi.org/project/XPRESSpipe
    :alt: PyPi Status

.. |doi| image:: https://zenodo.org/badge/######.svg
   :target: https://zenodo.org/badge/latestdoi/######
   :alt: DOI Status
