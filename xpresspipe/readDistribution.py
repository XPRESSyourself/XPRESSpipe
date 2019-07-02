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
    file1,
    file2=None):

    # Get files
    f1 = open(str(args_dict['input']) + str(file1), 'r').readlines()
    if file2 != None:
        f2 = open(str(args_dict['input']) + str(file1), 'r').readlines()

    # Get lengths of reads
    dist_list = Counter()
    for i, line in enumerate(f1[1:], 1):
        if (i - 1) % 4 == 0:
            length = len(line)
            if file2 != None:
                length += len(f2[i])
            dist_list[length] += 1

    # Clean up variables
    f1 = None
    if file2 != None:
        f2 = None

    # Compile length statistics
    distribution_profile = pd.DataFrame(dist_list, index=[0]).T
    distribution_profile.columns = ['count']

    # Export metrics
    distribution_profile['read size (bp)'] = distribution_profile.index
    distribution_profile.to_csv(
        str(args_dict['read_distributions']) + 'metrics/' + str(file1)[:-6] + '_metrics.txt',
        sep='\t')

"""Single-end RNA-seq pipeline"""
def se_dist(
        args):

    file, args_dict = args[0], args[1]

    output = str(file[8:-6]) # Get output file name before adding path to file name(s)

    get_distribution(args_dict, file)

"""Paired-end RNA-seq pipeline"""
def pe_dist(
        args):

    file1, file2, args_dict = args[0], args[1], args[2]

    # STAR first pass
    output = str(file1[8:-7]) # Get output file name before adding path to file name(s)

    get_distribution(args_dict, file1, file2)


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
    args_dict = add_directory(
        args_dict,
        'read_distributions',
        'individual_plots')

    # Get FASTQC file list and unzip
    files = get_files(
        args_dict['input'],
        ['.fastq', '.fq', '.txt'])

    if args_dict['type'] == 'PE':
        parallelize_pe(
            pe_dist,
            files,
            args_dict,
            mod_workers = True)
    else:
        parallelize(
            se_dist,
            files,
            args_dict,
            mod_workers = True)

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

    y = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_matrix_metrics(
            args_dict,
            str(args_dict['read_distributions']) + 'metrics/',
            file_list,
            'read size (bp)',
            'count',
            'read_distribution',
            args_dict['experiment'],
            args_dict['read_distributions'],
            str(args_dict['read_distributions']) + 'individual_plots/')

        y += 1
