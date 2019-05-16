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
from math import ceil

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_file_metrics
from .utils import add_directory, get_files
from .quality import run_fastqc

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
        'fastqc_output')
    args_dict = add_directory(
        args_dict,
        'read_distributions',
        'individual_plots')

    # Get FASTQC file list and unzip
    files = get_files(
        args_dict['input'],
        ['.fastq', '.fq', '.txt'])

    # Perform fastqc on each file and unzip output
    parallelize(
        run_fastqc,
        files,
        args_dict,
        mod_workers = True)

    files = get_files(
        args_dict['fastqc_output'],
        ['.zip'])

    for file in files:
        if file.endswith('.zip'):
            os.system(
                'unzip'
                + ' -n -q '
                + str(args_dict['fastqc_output']) + str(file)
                + ' -d ' + str(args_dict['fastqc_output'])
                + str(args_dict['log']))

    # Compile read distributions
    files = get_files(
        args_dict['fastqc_output'],
        ['.zip'])

    # Get metrics to plot
    files = [str(x[:-4]) + '/fastqc_data.txt' for x in files]

    file_number = ceil(len(files) / 6)
    file_lists = []

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+6])
        y += 6

    y = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_file_metrics(
            args_dict,
            args_dict['fastqc_output'],
            file_list,
            '#Length',
            '>>END_MODULE',
            'position',
            'read size (bp)',
            'fastqc',
            args_dict['experiment'],
            args_dict['read_distributions'],
            str(args_dict['read_distributions']) + 'individual_plots/')

        y += 1
