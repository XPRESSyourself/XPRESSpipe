###################
Building References
###################

Current curated reference files for compatible alignment programs are included
when installing RiboPipe for S. cerevisiae only due to file hosting size limits.
In order to use organisms other than yeast with RiboPipe, sub-modules are
included that will generate these curated reference files formatted for
compatibility with RiboPipe.


============================
Generating Curated Reference
============================

Currently, curated references compatible with RiboPipe can be generated for the
following model organisms:
Homo sapiens (STAR)

RiboPipe is easily adaptable to allow curation of other organisms as the need
arises. Future updates will allow for increased breadth of model organisms.

In order to generate these references, follow the example below:

1)  Create a directory where ever convenient that will be used to store references:

.. code-block:: shell

  $ cd ~/references

2)  Run RiboPipe's :data:`curate` sub-module, specifying references directory
created in previous step (or previously created for another curated reference):

.. code-block:: shell

  $ ribopipe curate -l /path/to/reference/ -p STAR -r human -c 16

If a model organism for a compatible alignment software is not available, please
contact us and we can aid in generating a different curated reference file.


=======================================
Using a Curated Reference with RiboPipe
=======================================

Once a curated reference has been generated, slightly different parameters must be
used with RiboPipe's :data:`riboseq` and :data:`rnaseq` sub-modules.

Provide the :data:`-c` or :data:`--custom` flags
Provide the full path to the reference folder with the :data:`-r` or :data:`--reference` flags,
such as in the example below:

.. code-block:: shell

  $ ribopipe riboseq -i /path/to/input -o /path/to/output ... \
  -c -r /scratch/general/lustre/u0690617/references/human_reference_star

If using alignment software other than STAR, the :data:`-p` or :data:`--program`
flags still need to be provided with the appropriate software name
