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
from functools import partial

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_file_metrics, compile_matrix_metrics, compile_complexity_metrics
from .utils import add_directory, get_files
from .gtfFlatten import flatten_reference, create_chromosome_index, create_coordinate_index
from .processBAM import read_bam, bam_sample
from .metaFeatures import meta_coordinates, psite_ranges, get_meta_profile, get_periodicity_profile

"""Get meta and periodicity indices from GTF"""
def get_indices(
        args_dict):

    # Read in GTF
    gtf = pd.read_csv(str(args_dict['reference']) + 'transcripts.gtf',
       sep='\t',
       header=None,
       comment='#',
       low_memory=False)

    # Flatten GTF
    gtf_flat = flatten_reference(
        gtf,
        threads = None)

    # Get GTF indices
    chromosome_index = create_chromosome_index(gtf_flat)
    coordinate_index = create_coordinate_index(gtf_flat)

    return chromosome_index, coordinate_index

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
        args,
        chromosome_index,
        coordinate_index):

    file, args_dict = args[0], args[1] # Parse args

    # Read in indexed bam file
    bam = read_bam(
        str(args_dict['input']) + str(file))

    bam_coordinates = psite_ranges(bam)
    if bam_coordinates == None:
        print('Warning: No reads passed filtering (28-30 nt) for file ' + str(file))
        return
  
    # Get profile
    profile_data = get_periodicity_profile(
        bam_coordinates,
        coordinate_index,
        chromosome_index)
    profile_data['index'] = profile_data.index
    profile_data.to_csv(
        str(args_dict['periodicity']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

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

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        ['_dedupRemoved.bam'])

    # Get indices
    chromosome_index, coordinate_index = get_indices(args_dict)

    # Perform periodicity analysis
    func = partial(
        get_peaks,
        chromosome_index = chromosome_index,
        coordinate_index = coordinate_index)
    parallelize(
        func,
        files,
        args_dict,
        mod_workers = True)

    # Get metrics to plot
    files = get_files(
        str(args_dict['periodicity']) + 'metrics/',
        ['_metrics.txt'])

    # Plot metrics for each file
    compile_matrix_metrics(
        args_dict,
        str(args_dict['periodicity']) + 'metrics/',
        files,
        'index',
        'phasing',
        'periodicity',
        args_dict['experiment'],
        args_dict['periodicity'])

"""Generate metagene profiles"""
def get_metagene(
        args,
        chromosome_index,
        coordinate_index):

    file, args_dict = args[0], args[1] # Parse args

    # Read in indexed bam file
    bam = read_bam(
        str(args_dict['input']) + str(file))
    bam_coordinates = meta_coordinates(bam)

    # Get profile
    profile_data = get_meta_profile(
        bam_coordinates,
        coordinate_index,
        chromosome_index)
    profile_data['index'] = profile_data.index
    profile_data.to_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

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

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        ['_dedupRemoved.bam'])

    # Get indices
    chromosome_index, coordinate_index = get_indices(args_dict)

    # Perform metagene analysis
    func = partial(
        get_metagene,
        chromosome_index = chromosome_index,
        coordinate_index = coordinate_index)
    parallelize(
        func,
        files,
        args_dict,
        mod_workers = True)

    # Compile metrics to plot
    files = get_files(
        str(args_dict['metagene']) + 'metrics/',
        ['_metrics.txt'])

    # Plot metrics for each file
    compile_matrix_metrics(
        args_dict,
        str(args_dict['metagene']) + 'metrics/',
        files,
        'index',
        'metacount',
        'metagene',
        args_dict['experiment'],
        args_dict['metagene'])

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
        args_dict['read_distributions'])

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
        ['dupRadar_metrics.txt'])

    # Plot metrics for each file
    compile_complexity_metrics(
        args_dict,
        str(args_dict['complexity']) + 'metrics/',
        files,
        'dupRateMulti',
        'RPKM',
        'library complexity (all_reads)',
        args_dict['experiment'],
        args_dict['complexity'])
