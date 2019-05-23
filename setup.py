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
__path__ = str(os.path.dirname(os.path.realpath(__file__))) + '/'

def test_system():

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
        develop.run(self)
        test_system()

class PostInstallCommand(install):
    # Post-installation for python setup.py install
    def run(self):
        install.run(self)
        test_system()

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
        'xpresstools',
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
