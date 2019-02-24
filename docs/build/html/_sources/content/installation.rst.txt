############
Installation
############

===================
Local Installation:
===================
1)  Make sure Python3, git, and wget are installed.
2)  Download `Conda <https://www.anaconda.com/download/#macos>`_, a package manager, for your operating system. Double click the `.pkg` file if on MacOS, the `.exe` file on Windows, or follow these `instructions <https://conda.io/docs/user-guide/install/linux.html#install-linux-silent>`_ on Linux.
3)  Execute the following lines of code in `Terminal <https://www.imore.com/how-use-terminal-mac-when-you-have-no-idea-where-start>`_ (on Mac, open Spotlight and type 'Terminal'):
4)  Download current repository:

.. code-block:: shell

  $ git clone https://github.com/j-berg/ribopipe.git
  $ cd ribopipe/ribopipe/references

5)  To download specific version

.. code-block:: shell

  $ tag='v0.1.4-beta'
  $ wget https://github.com/j-berg/ribopipe/archive/$tag.zip
  $ unzip ribopipe-${tag:1}.zip
  $ cd ribopipe-${tag:1}/ribopipe/references

6)  Get reference

.. code-block:: shell

  $ model='yeast'
  $ program='hisat2'
  $ wget https://sourceforge.net/projects/ribopipe/files/${program}_references/${model}_reference_${program}.zip
  $ unzip ${model}_reference_${program}.zip
  $ rm ${model}_reference_${program}.zip
  $ cd ../../
  $ python3 setup.py install --prefix ~/.local

7) add script installation location given near the end of the setup scripting output to ~/.bashrc or ~/.bash_profile
add to .bashrc

.. code-block:: shell

  $ echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bashrc

add to .bash_profile

.. code-block:: shell

  $ echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bash_profile

8) Test by typing the following:

.. code-block:: shell

  $ ribopipe --help

9) Install conda dependencies:

.. code-block:: shell

  $ ribopipe install

See local_install.sh in the `resources <https://github.com/j-berg/ribopipe/resources/>`_ folder for interactive script


=================
HPC Installation:
=================
1)  Make sure Python3, git, and wget are installed.
2)  Execute the following lines of code:

3)  To download current repository:

.. code-block:: shell

  $ git clone https://github.com/j-berg/ribopipe.git
  $ cd ribopipe/ribopipe/references

4)  To download specific version

.. code-block:: shell

  $ tag='v0.1.4-beta'
  $ wget https://github.com/j-berg/ribopipe/archive/$tag.zip
  $ unzip ribopipe-${tag:1}.zip
  $ cd ribopipe-${tag:1}/ribopipe/references

5)  Get reference

.. code-block:: shell

  $ model='yeast'
  $ program='hisat2'
  $ wget https://sourceforge.net/projects/ribopipe/files/${program}_references/${model}_reference_${program}.zip
  $ unzip ${model}_reference_${program}.zip
  $ rm ${model}_reference_${program}.zip
  $ cd ../../
  $ module load python3
  $ python setup.py install --prefix ~/.local

6)  Add script installation location given near the end of the setup scripting output to ~/.bashrc or ~/.bash_profile
add to .bashrc

.. code-block:: shell

  $ echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bashrc

add to .bash_profile

.. code-block:: shell

  $ echo "PATH='/path/to/scripts/:$PATH'" >> ~/.bash_profile

7) Test by typing the following:

.. code-block:: shell

  $ ribopipe --help

See hpc_install.sh in the `resources <https://github.com/j-berg/ribopipe/resources/>`_ folder for interactive script
