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

"""IMPORT DEPENDENCIES"""
import os
import sys
import pandas as pd
import numpy as np
import gc
from math import ceil
from functools import partial

"""IMPORT INTERNAL DEPENDENCIES"""
from .compile import compile_periodicity_metrics
from .processBAM import read_bam
from .utils import add_directory, get_files

transcriptome_alignment_read = 9
lower_quantile_bound = 0.125
upper_quantile_bound = 0.875

"""Manager for running periodicity summary plotting"""
def make_periodicity(
        args_dict):

    print('\nGenerating periodicity profiles...')
    args_dict = add_directory(
        args_dict,
        'output',
        'periodicity')
    args_dict = add_directory(
        args_dict,
        'periodicity',
        'metrics')

    # Get list of all files from input directory
    files = get_files(
        args_dict['input'],
        ['.bam'])


    # Get files read distributions
    for f in files:
        # Give non riboseq samples a bad suffix so they can't be searched by riboWaltz
        if f.endswith(args_dict['bam_suffix']):
            pass
        else:
            os.system(
                'mv'
                + ' ' + str(args_dict['input']) + str(f)
                + ' ' + str(args_dict['input']) + str(f) + '.rnaseq')

    if len(files) == 0:
        raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

    # Run riboWaltz in R
    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'Rperiodicity.r'
        + ' ' + str(args_dict['input'])
        + ' ' + str(args_dict['output'] + 'transcripts.idx')
        + ' ' + str(args_dict['periodicity']) + 'metrics/'
        + str(args_dict['log']))

    for f in files:
        if f.endswith(bam_suffix):
            pass
        else:
            os.system(
                'mv'
                + ' ' + str(args_dict['input']) + str(f) + '.rnaseq'
                + ' ' + str(args_dict['input']) + str(f))

    # Get metrics to plot
    files = get_files(
        str(args_dict['periodicity']) + 'metrics/',
        ['_metrics.txt'])

    file_number = ceil(len(files) / 6)
    file_lists = []

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+6])
        y += 6

    z = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_periodicity_metrics(
            str(args_dict['periodicity']) + 'metrics/',
            file_list,
            'periodicity_' + str(z),
            args_dict['experiment'],
            args_dict['periodicity'])

        z += 1

    chromosome_index = None
    coordinate_index = None
    gc.collect()
