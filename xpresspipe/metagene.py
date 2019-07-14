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
from math import ceil
import gc
from functools import partial

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_matrix_metrics
from .utils import add_directory, get_files
from .buildIndex import index_gtf

def roundup(value):
    return int(ceil(value))

def finish_metagene(args):

    file, args_dict = args[0], args[1]

    # Read in parsed bam
    bam = pd.read_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '.metaposit',
        sep = '\t')

    # Read in transcript exon space dictionary
    dict = pd.read_csv(
        str(args_dict['output']) + 'metagene.dict',
        sep = '\t')

    reference = pd.Series(dict['length'].values,index=dict['transcript']).to_dict()
    dict = None

    # Map exon space to each bam record based on its transcript ID
    bam['total_length']= bam['seqnames'].map(reference)
    bam = bam.dropna()

    # Calculate meta-position
    bam['meta_distance'] = bam['meta_position'] / bam['total_length'] * 100
    bam['meta_distance'] = bam['meta_distance'].apply(roundup)

    # Compile meta position statistics
    profile = pd.DataFrame()
    profile['metacount'] = [0 for x in range(1,101)]
    profile.index = range(1,101)

    for x in range(1,101):
        profile.loc[x] = bam[bam['meta_distance'] == x].shape[0]

    # Export metrics
    profile['representative transcript'] = profile.index
    profile.to_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '_metrics.txt',
        sep='\t')

    # Clean up
    bam = None
    reference = None
    dict = None
    profile = None
    gc.collect()
    os.system(
        'rm'
        + ' ' + str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '.metaposit')

def run_metagene(args):

    file, args_dict = args[0], args[1]
    file = '\"' + str(file) + '\"'

    print('Evaluating the metagene profile for ' + str(file))
    # Perform metagene analysis
    # Loop through each selected BAM
    # args[1] = Path to BAM files
    # args[2] = List of BAM files
    # args[3] = Index file with path
    # args[4] = Output file path
    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'Rmetagene.r'
        + ' ' + str(args_dict['input'])
        + ' ' + str(file)
        + ' ' + str(args_dict['metagene']) + 'metrics/'
        + str(args_dict['log']))

"""Manager for running metagene summary plotting"""
def make_metagene(
        args_dict):

    print('\nGenerating metagene profiles...')
    # Add output directory to output for metagene profiles
    args_dict = add_directory(
        args_dict,
        'output',
        'metagene')
    args_dict = add_directory(
        args_dict,
        'metagene',
        'metrics')

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        [str(args_dict['bam_suffix'])])

    # Get indices
    print('Generating index for genes...')
    index_gtf(
        args_dict,
        threads=None,
        geneCov=False,
        output=False)

    # Perform metagene analysis
    parallelize(
        run_metagene,
        files,
        args_dict,
        mod_workers = True)
    parallelize(
        finish_metagene,
        files,
        args_dict,
        mod_workers = True)
    os.system(
        'rm'
        + ' ' + str(args_dict['output']) + 'metagene.dict')

    # Compile metrics to plot
    print('Plotting...')
    files = get_files(
        str(args_dict['metagene']) + 'metrics/',
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
        compile_matrix_metrics(
            str(args_dict['metagene']) + 'metrics/',
            file_list,
            'representative transcript',
            'metacount',
            'metagene_' + str(z),
            args_dict['experiment'],
            args_dict['metagene'])
        z += 1

    gc.collect()
