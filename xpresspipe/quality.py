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
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_file_metrics, compile_matrix_metrics
from .utils import add_directory, get_files

"""Create MultiQC processing summary from all files in args_dict output"""
def get_multiqc_summary(
        args_dict):

    os.system(
        'multiqc'
        + ' ' + str(args_dict['output'])
        + ' -i ' + str(args_dict['experiment'])
        + ' -o ' + args_dict['output']
        + str(args_dict['log']))

"""Generate periodicity maps"""
def get_peaks(
        args):

    file, args_dict = args[0], args[1] # Parse args


"""Manager for running periodicity summary plotting"""
def make_periodicity(
        args_dict):

    args_dict = add_directory(
        args_dict,
        'output',
        'periodicity')
    args_dict = add_directory(
        args_dict,
        'periodicity',
        'metrics')
    output_location = args_dict['periodicity']

    # Get list of bam files from user input
    files = get_files(args_dict['input'], ['.bam'])

    # Perform metagene analysis
    parallelize(get_peaks, files, args_dict, mod_workers=True)

    # Get metrics to plot
    files = get_files(args_dict['metrics'], ['_periodicity_metagene_profile.txt'])

    # Plot metrics for each file
    compile_file_metrics(
        args_dict,
        str(output_location + 'metrics'),
        files,
        'metagene_average',
        str(int(args_dict['downstream'] - 1)),
        'position',
        'metagene coverage',
        'periodicity',
        args_dict['experiment'],
        output_location)

"""Generate metagene profiles"""
def get_profiles(
        args):

    file, args_dict = args[0], args[1] # Parse args

"""Manager for running metagene summary plotting"""
def make_metagene(
        args_dict):

    # Add output directory to output for metagene profiles
    args_dict = add_directory(
        args_dict,
        'output',
        'metagene')
    args_dict = add_directory(
        args_dict,
        'metagene',
        'metrics')
    output_location = args_dict['metagene']

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        ['.bam'])

    # Perform metagene analysis
    parallelize(
        get_profiles,
        files,
        args_dict,
        mod_workers = True)

    # Compile metrics to plot
    files = get_files(
        args_dict['metrics'],
        ['_rna_metrics'])

    # Plot metrics for each file
    compile_file_metrics(
        args_dict,
        str(output_location + 'metrics'),
        files,
        'normalized_position',
        None,
        'meta_position',
        'normalized_coverage (all_reads)',
        'metagene',
        args_dict['experiment'],
        output_location)

"""Generate read distribution profiles"""
def run_fastqc(
        args):

    file, args_dict = args[0], args[1]

    os.system(
        'fastqc'
        + ' -q ' + str(args_dict['input']) + str(file)
        + ' -o ' + str(args_dict['fastqc_output'])
        + str(args_dict['log']))

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
    output_location = args_dict['read_distributions']

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

    # Plot metrics for each file
    compile_file_metrics(
        args_dict,
        args_dict['fastqc_output'],
        files,
        '#Length',
        '>>END_MODULE',
        'position',
        'read size (bp)',
        'fastqc',
        args_dict['experiment'],
        output_location)

"""Measure library complexity"""
def run_complexity(
        args):

    file, args_dict = args[0], args[1]

    # Determine sequencing type
    if str(args_dict['type']).upper() == 'PE':
        paired = 'TRUE'
    else:
        paired = 'FALSE'

    # Run dupRadar in R
    os.system(
        'rscript'
        + ' ' + str(__path__) + '/complexity.r'
        + ' ' + str(args_dict['input']) + str(file)
        + ' ' + str(args_dict['gtf'])
        + ' ' + str(paired)
        + ' ' + str(args_dict['threads'])
        + ' ' + str(args_dict['metrics']) + str(file[:-4]) + 'dupRadar_metrics.txt'
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
    output_location = args_dict['complexity']

    # Get BAM files
    files = get_files(
        args_dict['input'],
        ['_deduplicated.bam'])

    # Perform metagene analysis
    parallelize(
        run_complexity,
        files,
        args_dict,
        mod_workers = False)

    # Get metrics to plot
    files = get_files(
        args_dict['metrics'],
        ['dupRadar_metrics.txt'])

    # Plot metrics for each file
    compile_matrix_metrics(
        args_dict,
        str(output_location + 'metrics'),
        files,
        'dupRateMulti',
        'RPKM',
        'library complexity (all_reads)',
        args_dict['experiment'],
        output_location)
