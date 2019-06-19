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
import gc
from math import ceil
from functools import partial

"""IMPORT INTERNAL DEPENDENCIES"""
from .parallel import parallelize
from .compile import compile_matrix_metrics
from .utils import add_directory, get_files
from .processBAM import read_bam, bam_sample
from .quality import get_indices, get_position

"""P site ranges for 28-30mers"""
def psite_ranges(
        bam):

    # Keep only optimal footprint size
    bam = bam[(bam[9].str.len() == 28)]

    # Get rightmost coordinates for each read
    bam[16] = bam[3] + bam[9].str.len()

    # Return as array
    phased_coordinates = bam[[2,3,16]]

    if len(phased_coordinates.index) == 0:
        return None
    else:
        return phased_coordinates.values.tolist()

"""Search for plausible records per provided coordinate"""
def get_coordinate_records_period(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate_reverse,
        search_coordinate_forward):

    search_coordinate_reverse += 16
    search_coordinate_forward -= 16

    record_array = []
    chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
    for index, record in enumerate(chromosome_array):
        if record[0] <= search_coordinate_reverse and record[1] >= search_coordinate_reverse \
        or record[0] <= search_coordinate_forward and record[1] >= search_coordinate_forward:
            record_array.append(record)

    return record_array

"""Get periodicity profile for bam file"""
def get_periodicity_profile(
        aligned_reads_index,
        coordinate_index,
        chromosome_index,
        start_range=-1,
        end_range=101):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(
        0,
        index = range(start_range, end_range),
        columns = ['count'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records_period(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1],
                record[2])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            position_count = []
            for index, transcript_record in enumerate(record_array):
                # Determine P-site coordinate to search
                if transcript_record[1] == '+':
                    coordinate = record[2] - 16
                else:
                    coordinate = record[1] + 16

                # Get position relative to start of the p-site
                position = get_position(
                    coordinate,
                    transcript_record[3],
                    transcript_record[2])

                if position != None and position > start_range and position < end_range:
                    # Position cannot be inclusive of start or end coorindate or will throw invalid index error
                    position_count.append(position)

            for x in position_count:
                count = 1 / len(position_count)
                metagene_profile.at[x, 'count'] += count

    return metagene_profile

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
        print('Warning: No reads passed filtering (28 nt) for file ' + str(file))
        return

    bam = None # Some clean-up
    del bam
    gc.collect()

    # Get profile
    profile_data = get_periodicity_profile(
        bam_coordinates,
        coordinate_index,
        chromosome_index)
    profile_data['position from start'] = profile_data.index
    profile_data.to_csv(
        str(args_dict['periodicity']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

"""Manager for running periodicity summary plotting"""
def make_periodicity(
        args_dict):

    print('\nGenerating periodicity profiles...')
    args_dict = add_directory(
        args_dict,
        'output',
        'periodicity')
    args_dict = add_directory(
        args_dict,
        'periodicity',
        'metrics')
    args_dict = add_directory(
        args_dict,
        'periodicity',
        'individual_plots')

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        [str(args_dict['bam_suffix'])])

    # Get indices
    chromosome_index, coordinate_index = get_indices(args_dict, record_type='CDS')

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
            str(args_dict['periodicity']) + 'metrics/',
            file_list,
            'position from start',
            'count',
            'periodicity',
            args_dict['experiment'],
            args_dict['periodicity'],
            str(args_dict['periodicity']) + 'individual_plots/')

        y += 1

    chromosome_index = None
    coordinate_index = None
    gc.collect()
