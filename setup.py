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

from setuptools import setup
import re

with open('xpresspipe/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

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
    install_requires=[
          
      ],

    entry_points = {
        'console_scripts': [
            'xpresspipe = xpresspipe.__main__:main'
        ]
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)
