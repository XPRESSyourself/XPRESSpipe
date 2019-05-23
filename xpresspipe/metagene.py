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
from .processBAM import read_bam, bam_sample
from .quality import get_indices, get_position

"""Retrieve meta-coordinate for each read in input BAM matrix"""
def meta_coordinates(
        bam):

    # Map middle point of each read as left-most position plus half of read length
    bam[16] = bam[3] + (bam[9].str.len() / 2).apply(np.floor).astype('int64')

    # Return as array
    mid_coordinates = bam[[2,16]]
    return mid_coordinates.values.tolist()

"""Search for plausible records per provided coordinate"""
def get_coordinate_records_meta(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate_start):

    record_array = []
    chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
    for index, record in enumerate(chromosome_array):
        if record[0] <= search_coordinate_start and record[1] >= search_coordinate_start:
            record_array.append(record)

    return record_array

"""Get meta profile for bam file"""
def get_meta_profile(
        aligned_reads_index,
        coordinate_index,
        chromosome_index):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(
        0,
        index = range(101),
        columns = ['metacount'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records_meta(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            position_count = []
            for index, transcript_record in enumerate(record_array):
                position = get_position(
                    record[1],
                    transcript_record[3],
                    transcript_record[2])
                if position != None:
                    position_count.append([position, transcript_record[4]])

            for x in position_count:
                meta_position = int((x[0] / x[1]) * 100)
                count = 1 / len(position_count)
                metagene_profile.at[meta_position, 'metacount'] += count

    return metagene_profile

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
    profile_data['meta-transcript'] = profile_data.index
    profile_data.to_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

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
    args_dict = add_directory(
        args_dict,
        'metagene',
        'individual_plots')

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

    file_number = ceil(len(files) / 6)
    file_lists = []

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+6])
        y += 6

    y = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_matrix_metrics(
            args_dict,
            str(args_dict['metagene']) + 'metrics/',
            file_list,
            'meta-transcript',
            'metacount',
            'metagene',
            args_dict['experiment'],
            args_dict['metagene'],
            str(args_dict['metagene']) + 'individual_plots/')

        y += 1

    chromosome_index = None
    coordinate_index = None
    gc.collect()
