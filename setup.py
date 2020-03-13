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
import re
import os
import sys
import subprocess
from setuptools import setup
__path__ = str(os.path.dirname(os.path.realpath(__file__))) + '/'

def get_cufflinks():

    # Get python version
    python =  sys.version_info[0]

    if 'darwin' in sys.platform:
        system = 'MacOSX'
        cufflinks = 'cufflinks-2.1.1.OSX_x86_64'

    elif 'linux' in sys.platform:
        system = 'Linux'
        cufflinks = 'cufflinks-2.1.1.Linux_x86_64'

    else:
        raise Exception('Cannot recognize operating system. Expected \"darwin\" or \"linux\", got ' + str(sys.platform))

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

def move_fastp_lite():

    # Check fastp_lite exec is there first
    

    subprocess.call(
        'mv ' + str(__path__) + 'fastp_lite/fastp_lite ' + str(__path__) + 'xpresspipe; '
        + 'echo "fastp_lite installed"; '),
        shell = True)

"""Get version"""
with open('xpresspipe/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

"""Setup arguments"""
get_cufflinks()
move_fastp_lite()
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

    entry_points = {
        'console_scripts': [
            'xpresspipe = xpresspipe.__main__:main'
        ]
    },

    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)
