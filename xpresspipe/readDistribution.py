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
from math import ceil
from collections import Counter

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize, parallelize_pe
from .compile import compile_matrix_metrics
from .utils import add_directory, get_files


"""Get distributions matrix for SE or PE files"""
def get_distribution(
    args_dict,
    output,
    file1,
    file2=None):

    print('Evaluating the read distribution profile for ' + str(file1))
    # Perform metagene analysis
    # Loop through each selected BAM
    # args[1] = Path to BAM files
    # args[2] = List of BAM files
    # args[3] = Index file with path
    # args[4] = Output file path
    os.system(
        'julia'
        + ' ' + str(args_dict['path']) + 'readDistribution.jl'
        + ' ' + str(file1)
        + ' ' + str(file2)
        + ' ' + str(output)
        + str(args_dict['log']))

"""Format output file name
"""
def format_name(
        file):

    if file[-3:] == '.fq':
        file = file[:-3]

    elif file[-6:] == '.fastq':
        file = file[:-6]

    else:
        pass

    if file[:8] == 'trimmed_':
        file = file[8:]

    else:
        pass

    return file

"""Single-end RNA-seq pipeline"""
def se_dist(
        args):

    file, args_dict = args[0], args[1]

    file_input = str(args_dict['input']) + str(file)

    # Get output file name before adding path to file name(s)
    file_output = format_name(file)
    output = str(args_dict['read_distributions']) + 'metrics/' + str(file_output) + '_metrics.txt'

    get_distribution(args_dict, output, file_input)

"""Paired-end RNA-seq pipeline"""
def pe_dist(
        args):

    file1, file2, args_dict = args[0], args[1], args[2]

    file1_input = str(args_dict['input']) + str(file1)
    file2_input = str(args_dict['input']) + str(file2)

    # Get output file name before adding path to file name(s)
    file_output = format_name(file1)
    output = str(args_dict['read_distributions']) + 'metrics/' + str(file_output) + '_metrics.txt'

    get_distribution(args_dict, output, file1_input, file2_input)


"""Manager for running read distribution summary plotting"""
def make_readDistributions(
        args_dict):

    args_dict = add_directory(
        args_dict,
        'output',
        'read_distributions')
    args_dict = add_directory(
        args_dict,
        'read_distributions',
        'metrics')

    # Get FASTQC file list and unzip
    files = get_files(
        args_dict['input'],
        ['.fastq', '.fq', '.txt'])

    if args_dict['type'] == 'PE':
        parallelize_pe(
            pe_dist,
            files,
            args_dict,
            mod_workers = 'all')
    else:
        parallelize(
            se_dist,
            files,
            args_dict,
            mod_workers = 'all')

    # Get metrics to plot
    files = get_files(
        str(args_dict['read_distributions']) + 'metrics/',
        ['_metrics.txt'])

    file_number = ceil(len(files) / 6)
    file_lists = []

    if file_number < 1:
        file_number = 1

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+6])
        y += 6

    z = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_matrix_metrics(
            str(args_dict['read_distributions']) + 'metrics/',
            file_list,
            'read size (bp)',
            'count',
            'read_distribution_' + str(z),
            args_dict['experiment'],
            args_dict['read_distributions'])

        z += 1
