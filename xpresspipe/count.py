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
import os
import sys
import datetime
import pandas as pd
from xpresstools import count_table

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, add_directory
from .parallel import parallelize

"""Compile counts tables from HTseq output files"""
def count_file(args):

    file, args_dict = args[0], args[1] # Parse args

    # Count
    os.system('htseq-count'
        + ' -q'
        + ' -m intersection-nonempty'
        + ' -t exon'
        + ' -i gene_id'
        + ' -r pos'
        + ' -s no'
        + ' ' + str(args_dict['input']) + str(file)
        + ' ' + str(args_dict['gtf'])
        + ' > ' + str(args_dict['counts']) + str(file[:-4]) + '.tsv')

"""Run count reads manager"""
def count_reads(args_dict):

    # Add output directories
    args_dict = add_directory(args_dict, 'output', 'counts')

    # Get list of files to count based on acceptable file types
    files = get_files(args_dict['input'], ['.bam'])

    # Count aligned RNAseq reads
    parallelize(count_file, files, args_dict, mod_workers=True)

    return args_dict

"""Take directory of single counts files and collate into single table"""
def collect_counts(args_dict):

    # Add output directories
    args_dict = add_directory(args_dict, 'output', 'counts')

    # Get list of files to count based on acceptable file types
    files = get_files(args_dict['input'], ['.tsv'])

    # Append path to file list
    count_files = []
    for x in files:
        count_files.append(str(args_dict['input']) + str(x))

    # Create and output collated count table
    counts = count_table(count_files)

    # Output data
    if 'experiment' in args_dict and args_dict['experiment'] != None:
        counts.to_csv(str(args_dict['counts']) + str(args_dict['experiment']) + '_counts_table.tsv', sep='\t')
    else:
        cdt = datetime.datetime.now()
        counts.to_csv(str(args_dict['counts']) + str(cdt.year) + '_' + str(cdt.month) + '_' + str(cdt.day) + '_' + str(cdt.hour) + 'h_' + str(cdt.minute) + 'm_' + str(cdt.second) + 's_counts_table.tsv', sep='\t')
