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

"""IMPORT INTERNAL DEPENDENCIES"""
from .gtfFlatten import flatten_reference, create_chromosome_index, create_coordinate_index, make_flatten
from .utils import add_directory, get_files
from .parallel import parallelize

"""Get meta and periodicity indices from GTF"""
def get_indices(
    args_dict,
    record_type='exon',
    gene_name=None,
    threads=None):

    # Read in GTF
    gtf = pd.read_csv(
        str(args_dict['gtf']),
        sep='\t',
        header=None,
        comment='#',
        low_memory=False)

    if gene_name != None:
        gtf = gtf.loc[gtf[8].str.contains(str(gene_name))]
        gtf = gtf.reset_index(drop=True)

    # Flatten GTF
    if args_dict['gtf'].endswith('_LC.gtf') == True:
        gtf_flat = make_flatten(
            gtf,
            record_type)
    else:
        gtf_flat = flatten_reference(
            gtf,
            record_type,
            threads = threads)

    # Get GTF indices
    chromosome_index = create_chromosome_index(gtf_flat)
    coordinate_index = create_coordinate_index(gtf_flat)

    return chromosome_index, coordinate_index

"""Get relative position of coordinate to the start of the transcript"""
def get_position(
        position,
        coordinates,
        strand):

    # Can order operations same since reference puts - strand records in reverse already
    location = 0
    last_coordinate = 0

    for y in coordinates:

        # Exit if mapped to intron
        if strand == '+':
            next_coordinate = min(y)
        else: # '-'
            next_coordinate =  max(y)

        if position > last_coordinate and position < next_coordinate:
            return None

        # Map to exon position
        if position >= min(y) and position <= max(y):
            if strand == '+':
                location += abs(position - min(y))
            else: # '-'
                location += abs(max(y)- position)

            return location

        else:
            location += abs(y[1] - y[0])
            if strand == '+':
                next_coordinate = max(y)
            else: # '-'
                next_coordinate =  min(y)

"""Generate read distribution profiles"""
def run_fastqc(
        args):

    file, args_dict = args[0], args[1]

    os.system(
        'fastqc'
        + ' -q ' + str(args_dict['input']) + str(file)
        + ' -o ' + str(args_dict['fastqc'])
        + str(args_dict['log']))

"""Parellel run of FastQC"""
def get_fastqc(args_dict):

    args_dict = add_directory(
        args_dict,
        'output',
        'fastqc')

    files = get_files(
        args_dict['input'],
        ['.fastq', '.fq', '.txt'])

    # Perform fastqc on each file and unzip output
    parallelize(
        run_fastqc,
        files,
        args_dict,
        mod_workers = True)

"""Create MultiQC processing summary from all files in args_dict output"""
def get_multiqc_summary(
        args_dict):

    os.system(
        'multiqc'
        + ' ' + str(args_dict['output'])
        + ' -i ' + str(args_dict['experiment'])
        + ' -o ' + args_dict['output']
        + str(args_dict['log']))
