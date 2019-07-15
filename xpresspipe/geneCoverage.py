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
from .compile import compile_coverage
from .utils import add_directory, get_files
from .buildIndex import index_gtf

plots_per_page = 8

def run_coverage(args):

    file, args_dict = args[0], args[1]
    file = '\"' + str(file) + '\"'

    print('Evaluating the gene coverage of ' + str(file))
    # Perform metagene analysis
    # Loop through each selected BAM
    # args[1] = Path to BAM files
    # args[2] = List of BAM files
    # args[3] = Index file with path
    # args[4] = Output file path
    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'RgeneCoverage.r'
        + ' ' + str(args_dict['input'])
        + ' ' + str(file)
        + ' ' + str(args_dict['output']) + str(args_dict['gene_name']) + '.idx'
        + ' ' + str(args_dict['coverage']) + 'metrics/' + str(args_dict['gene_name']) + '_'
        + str(args_dict['log']))

"""
Func: Get coverage profile for a specific gene
- Creates output directories for coverage (parent) and metrics
- Gets BAM files to plot coverage profiles
- Optionally orders files if specified by user, if not will order alphanumerically
- Generates chromosome and coordinate indices for gene of interest
- Generates coverage tables for each file
- Plots each file as a summary
@param args_dict: Global user arguments dictionary
@return: None, pipes to plotting of output metrics
"""
def make_coverage(
    args_dict):

    # Add output directories
    print('\nGenerating gene coverage profiles...')
    args_dict = add_directory(
        args_dict,
        'output',
        'coverage')

    args_dict = add_directory(
        args_dict,
        'coverage',
        'metrics')

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        [str(args_dict['bam_suffix'])])
    if len(files) == 0:
        raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

    # Get samples user specified
    if args_dict['samples'] != None:
        sample_list = []
        for x in args_dict['samples']:
            for y in files:
                if y in x:
                    sample_list.append(y)
                    break
        files = sample_list

    # Perform gene coverage analysis
    parallelize(
        run_coverage,
        files,
        args_dict,
        mod_workers = True)

    # Compile metrics to plot
    print('Plotting...')
    regions = pd.read_csv(
        str(args_dict['output']) + str(args_dict['gene_name']) + '.fts',
        sep='\t')

    files = get_files(
        str(args_dict['coverage']) + 'metrics/',
        ['_metrics.txt'])
    files = [f for f in files if f.startswith(str(args_dict['gene_name']))]

    page_number = ceil(len(files) / plots_per_page)
    file_lists = []

    y = 0
    for x in range(page_number):
        file_lists.append(files[y:y+plots_per_page])
        y += plots_per_page

    z = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_coverage(
            str(args_dict['coverage']) + 'metrics/',
            file_list,
            args_dict['gene_name'],
            regions,
            args_dict['sample_names'],
            str(args_dict['gene_name']) + '_geneCoverage' + str(z),
            args_dict['coverage'],
            args_dict['plot_color'])

        z += 1

    gc.collect()
