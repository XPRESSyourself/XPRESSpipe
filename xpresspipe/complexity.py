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

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_complexity_metrics
from .utils import add_directory, get_files

"""Measure library complexity"""
def run_complexity(
        args):

    file, args_dict = args[0], args[1]

    # Determine sequencing type
    if str(args_dict['type']).upper() == 'PE':
        paired = 'True'
    else:
        paired = 'False'

    # Run dupRadar in R
    os.system(
        'rscript'
        + ' ' + str(args_dict['path']) + '/Rcomplexity.r'
        + ' ' + str(args_dict['input']) + str(file)
        + ' ' + str(args_dict['gtf'])
        + ' ' + str(paired)
        + ' ' + str(args_dict['threads'])
        + ' ' + str(args_dict['complexity']) + 'metrics/' + str(file[:-4]) + '_metrics.txt'
        + str(args_dict['log']))

"""Manager for running complexity summary plotting"""
def make_complexity(args_dict):

    args_dict = add_directory(
        args_dict,
        'output',
        'complexity')
    args_dict = add_directory(
        args_dict,
        'complexity',
        'metrics')

    # Get BAM files
    files = get_files(
        args_dict['input'],
        ['_dedupMarked.bam'])

    # Perform metagene analysis
    parallelize(
        run_complexity,
        files,
        args_dict,
        mod_workers = False)

    # Get metrics to plot
    files = get_files(
        str(args_dict['complexity']) + 'metrics/',
        ['_metrics.txt'])

    # Plot metrics for each file
    compile_complexity_metrics(
        args_dict,
        str(args_dict['complexity']) + 'metrics/',
        files,
        'RPK',
        'dupRate',
        'library_complexity',
        args_dict['experiment'],
        args_dict['complexity'])
