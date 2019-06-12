"""
XPRESSpipe
An alignment and analysis pipeline for RNAseq data
alias: xpresspipe

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

"""IMPORT DEPENDENCIES"""
from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install
import re
import os
import sys
import subprocess

__path__ = str(os.path.dirname(os.path.realpath(__file__))) + '/'

"""Test system for cufflinks compatibility"""
def test_system():

    # Get python version
    python =  sys.version_info[0]

    if 'darwin' in sys.platform:
        system = 'MacOSX'
        cufflinks = 'cufflinks-2.1.1.OSX_x86_64'

    elif 'linux' in sys.platform:
        system = 'Linux'
        cufflinks = 'cufflinks-2.1.1.Linux_x86_64'

    else:
        raise Exception('Cannot recognize operating system')

    # Install Cufflinks
    subprocess.call(
        ('echo "Installing cufflinks binary for ' + str(system) + '..."; '
        + 'curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/' + str(cufflinks) + '.tar.gz -o ' + str(__path__) + 'cufflinks.tar.gz; '
        + 'tar -zxvf ' + str(__path__) + 'cufflinks.tar.gz; '
        + 'rm ' + str(__path__) + 'cufflinks.tar.gz; '
        + 'mv ' + str(__path__) + str(cufflinks) + ' ' + str(__path__) + 'cufflinks; '
        + 'mv ' + str(__path__) + 'cufflinks/cufflinks ' + str(__path__) + 'xpresspipe; '
        + 'echo "Cufflinks installed"; '),
        shell = True)

    # Install Conda
    if 'conda' not in sys.version:
        conda_input = input(
            'It does not appear the package manager, Conda, is installed on your system.'
            + ' XPRESSpipe uses this package manager to install required dependencies of this software.'
            + ' If you would like to install this manually, please visit https://conda.io/projects/conda/en/latest/user-guide/install/index.html.'
            + '\nWould you like XPRESSpipe to download this software for you? [y/n]: ')
        if conda_input.lower() == 'y' or conda_input.lower() == 'yes':
            subprocess.call(
                ('echo "Installing conda for ' + str(system) + '..."; '
                + 'curl https://repo.anaconda.com/miniconda/Miniconda' + str(python) + '-latest-' + str(system) + '-x86_64.sh -o ~/miniconda.sh; '
                + 'bash ~/miniconda.sh -b -p ~/miniconda; '
                + 'conda config --set always_yes yes --set changeps1 no --set show_channel_urls true; '
                + 'conda config --add channels defaults; '
                + 'conda config --add channels r; '
                + 'conda config --add channels bioconda; '
                + 'conda config --add channels conda-forge; '
                + 'echo "Conda installed"; '),
                shell = True)

    # Install Conda Dependencies
    install_input = input(
        'Install conda dependencies via XPRESSpipe setup? [y/n]: ')

    if install_input.lower() == 'y' or install_input.lower() == 'yes':
        subprocess.call(
            ('echo "Installing conda dependencies..."; '
            + 'conda config --add channels defaults; '
            + 'conda config --add channels r; '
            + 'conda config --add channels bioconda; '
            + 'conda config --add channels conda-forge; '
            + 'conda install -y fastp STAR samtools bedtools deeptools fastqc htseq pandas numpy biopython scipy r conda-forge::ncurses libiconv bioconductor-rsubread bioconductor-dupradar bioconductor-deseq2 matplotlib=2.2.3; '
            + 'echo "Conda dependencies installed"; '),
            shell = True)

"""Get version"""
with open('xpresspipe/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

"""Setup arguments"""
test_system()
setup(
    name = 'XPRESSpipe',
    version = version,
    description = 'An alignment and analysis pipeline for RNAseq data',
    long_description = open('README.md').read(),
    author = 'Jordan Berg',
    author_email = 'jordan.berg@biochem.utah.edu',
    url = 'https://github.com/XPRESSyourself/XPRESSpipe',
    packages = ['xpresspipe'],
    exclude = ['tests','docs'],
    package_dir = {'xpresspipe': 'xpresspipe'},
    include_package_data = True,
    license = 'GPL-3.0',
    zip_safe = False,

    install_requires = [
        'xpressplot',
        'multiqc',
        'numpydoc',
        'psutil'
    ],

    entry_points = {
        'console_scripts': [
            'xpresspipe = xpresspipe.__main__:main'
        ]
    },

    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)
