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

import sys

def build_command():
    print('The following program will help you design a command for use with XPRESSpipe\n')
    pipeline = -1
    while pipeline == -1:

        option = input('What do you want to do:\n1: Curate reference\n2: Process data\n3: Analyze data\n')
        if option == '1':
            print('coming soon...')
        elif option == '2':
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
            indir = input('Where is your raw sequence data found? (provide full path to directory with .fastq files)\n')

            # Get output directory
            outdir = input('Where would you like your processed data output to? (provide full path)\n')

            # Get adaptors
            adaptors_bool = input('Were adaptors used in generating your sequence libraries? (yes/no)\n')
            if adaptors_bool.lower() == 'yes':
                polyx_bool = input('Were polyX adaptors used? (yes/no)\n')
                if polyx_bool.lower() == 'yes':
                    if pipeline == 'seRNAseq' or pipeline == 'riboseq':
                        adaptors = "polyX"
                    else:
                        raise Warning('polyX adaptors not compatible with paired-end sequencing module currently. Please try again')
                        pipeline = -1
                else:
                    if pipeline == 'seRNAseq' or pipeline == 'riboseq':
                        adaptors = input('Please provide the adaptor used:\n')
                    else:
                        adaptors = input('Please provide the adaptors used separated by a space:\n')
                        if ' ' not in adaptors:
                            raise Warning('Input did not include a space. Please try again')
                            pipeline = -1
            else:
                if pipeline == 'peRNAseq':
                    adaptors = 'None None'
                else:
                    adaptors = "None"


            # Summarize
            print('\nCommand summary:')
            print('Pipeline = ' + str(pipeline))
            print('Input Path = ' + str(indir))
            print('Output Path = ' + str(outdir))

            print('\nHere is the command you should use to process your data:\n')
            print('xpresspipe ' + str(pipeline) + ' -i ' + str(indir) + ' -o ' + str(outdir) + '\n')
            sys.exit(1)
        elif option == '3':
            print('coming soon...')
        else:
            print('Invalid input, try again...\n')
            pipeline = -1
