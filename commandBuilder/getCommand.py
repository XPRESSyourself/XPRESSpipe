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
from __future__ import print_function
from flask import Flask

app = Flask(__name__)



print('The following program will help you design a command for use with XPRESSpipe\n')

# Get pipeline
pipeline = -1
while pipeline == -1:

    pipeline = input('What kind of processing would you like to perform (select the number)?\n1: Single-end RNA-seq\n2: Paired-end RNA-seq\n3: Ribosome Profiling\n')

    if pipeline == '1':
        pipeline = 'seRNAseq'
    elif pipeline == '2':
        pipeline = 'peRNAseq'
    elif pipeline == '3':
        pipeline = 'riboseq'
    else:
        print('Invalid input, try again...\n')
        pipeline = -1

# Get input directory
indir = input('Specify input directory (provide full path)\n')

# Get output directory
outdir = input('Specify output directory (provide full path)\n')

















# Summarize
print('\nCommand summary:')
print('Pipeline = ' + str(pipeline))
print('Input Path = ' + str(indir))
print('Output Path = ' + str(outdir))

print('\nHere is the command you should use to process your data:')
print('xpresspipe ' + str(pipeline) + ' -i ' + str(indir) + ' -o ' + str(outdir))
