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

"""Test system for cufflinks compatibility"""

def get_path():

    arguments = sys.argv

    if arguments.count('--prefix') > 1 or arguments.count('--install_dir') > 1:
        raise Exception('Argument inputs are incorrect. You provided too many instances of either --prefix or --install_dir')

    for x in range(len(arguments)):
        if arguments[x] == '--prefix' or arguments[x] == '--install_dir':
            return arguments[x + 1]

def test_system(__path__):

    if sys.platform == 'darwin':
        subprocess.call(
            ('echo "Installing cufflinks binary for macOS..."; '
            + 'curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.OSX_x86_64.tar.gz -o ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'tar -zxvf ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'rm ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'mv ' + str(__path__) + 'cufflinks-2.1.1.OSX_x86_64 ' + str(__path__) + 'cufflinks; '
            + 'echo "Cufflinks installed"; '),
            shell = True)

    elif sys.platform == 'linux':
        subprocess.call(
            ('echo "Installing cufflinks binary for Linux..."; '
            + 'curl http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz -o ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'tar -zxvf ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'rm ' + str(__path__) + 'cufflinks.tar.gz; '
            + 'mv ' + str(__path__) + 'cufflinks-2.1.1.Linux_x86_64 ' + str(__path__) + 'cufflinks; '
            + 'echo "Cufflinks installed"; '),
            shell = True)

    else:
        raise Exception('Cannot recognize operating system')

"""Define post-install classes"""
class PostDevelopCommand(develop):
    # Post-installation for python setup.py develop
    def run(self):
        __path__ = get_path()
        develop.run(self)
        __path__ = get_path()
        test_system(__path__)

class PostInstallCommand(install):
    # Post-installation for python setup.py install
    def run(self):
        install.run(self)
        __path__ = get_path()
        test_system(__path__)

"""Get version"""
with open('xpresspipe/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

"""Setup arguments"""
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
        'numpydoc'

    ],

    cmdclass = {
        'develop': PostDevelopCommand,
        'install': PostInstallCommand
    },

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
