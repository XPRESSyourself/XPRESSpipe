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
from .compile import compile_matrix_metrics, compile_coverage
from .utils import add_directory, get_files
from .processBAM import read_bam, bam_sample
from .quality import get_indices, get_position

"""
Func: Retrieve meta-coordinate for each read in input BAM matrix
- Find middle coordinate of each read
@param bam: BAM file in pandas dataframe format
@return: Middle coordinates and respective chromosome as multi-dimensional list
"""
def meta_coordinates(
        bam):

    # Map middle point of each read as left-most position plus half of read length
    bam[16] = bam[3] + (bam[9].str.len() / 2).apply(np.floor).astype('int64')

    # Return as array
    mid_coordinates = bam[[2,16]]
    return mid_coordinates.values.tolist()

"""
Func: Search for plausible records per provided coordinate
- Searches for a gene range where search coordinate falls within
@param coordinate_index: Multi-dimensional array for each gene with start, stop, strand, length, and exon region information
@param chromosome_index: Dictionary with chromosome number and ordered position in the index
@param search_chromosome: Chromosome number of coordinate to search
@param search_coordinate_start: Leftmost coordinate of search read
@return: Updated record_array if search coordinate falls in gene range
"""
def get_coordinate_records_meta(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate_start):

    record_array = []
    try:
        chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
        for index, record in enumerate(chromosome_array):
            if record[0] <= search_coordinate_start and record[1] >= search_coordinate_start:
                record_array.append(record)

        return record_array
    except:
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

    bam = None # Some clean-up
    del bam

    # Get profile
    profile_data = get_meta_profile(
        bam_coordinates,
        coordinate_index,
        chromosome_index)
    profile_data['representative transcript'] = profile_data.index

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
        [str(args_dict['bam_suffix'])])

    # Get indices
    chromosome_index, coordinate_index = get_indices(args_dict, record_type=args_dict['type'])

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

    for file_list in file_lists:

        # Plot metrics for each file
        compile_matrix_metrics(
            args_dict,
            str(args_dict['metagene']) + 'metrics/',
            file_list,
            'representative transcript',
            'metacount',
            'metagene',
            args_dict['experiment'],
            args_dict['metagene'],
            str(args_dict['metagene']) + 'individual_plots/')

    chromosome_index = None
    coordinate_index = None
    gc.collect()















def cov_coordinates(
    bam):

    # Map middle point of each read as left-most position plus half of read length
    bam[16] = bam[3] + (bam[9].str.len()).astype('int64')

    # Return as array
    mid_coordinates = bam[[2,3,16]]
    return mid_coordinates.values.tolist()


def get_cov_position(
    position,
    coordinates):

    # Can order operations same since reference puts - strand records in reverse already
    location = 0
    last_coordinate = 0

    for y in coordinates:

        # Map to exon position
        if position >= min(y) and position <= max(y):
            return position
        else:
            pass

    return None

"""Get meta profile for bam file"""
def get_cov_profile(
    aligned_reads_index,
    coordinate_index,
    chromosome_index):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(
        0,
        index = range(coordinate_index[0][0][0], coordinate_index[0][0][1] + 1),
        columns = ['raw_count'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        for x in range(record[1], (record[2] + 1)):
            position = get_cov_position(
                x,
                coordinate_index[0][0][3])

            if position != None:
                metagene_profile.at[position, 'raw_count'] += 1

    return metagene_profile

"""
func: Generate metagene profiles
- Parse through BAM file for plausible reads falling within specified gene
"""
def get_coverage(
    args,
    chromosome_index,
    coordinate_index):

    file, args_dict = args[0], args[1] # Parse args

    # Read in indexed bam file
    bam = read_bam(
        str(args_dict['input']) + str(file))

    # Get BAM file relevant to gene of interest if using geneCoverage
    chr = []
    for key, value in chromosome_index.items():
        chr.append(key)

    if len(coordinate_index) == 1 and len(coordinate_index[0]) == 1 and len(chr) == 1:
        bam = bam.loc[(bam[2] == chr[0]) & (bam[3] >= (coordinate_index[0][0][0] - coordinate_index[0][0][4])) & (bam[3] <= (coordinate_index[0][0][1] + coordinate_index[0][0][4]))]
    else:
        raise Exception('A single gene record was not selected.')

    # Get meta-data for bam file relevant to gene region
    bam_coordinates = cov_coordinates(bam)

    bam = None # Some clean-up
    del bam

    # Get profile
    profile_data = get_cov_profile(
        bam_coordinates,
        coordinate_index,
        chromosome_index)
    profile_data['transcript'] = profile_data.index

    coding_space = []
    for x in coordinate_index[0][0][3]:
        for y in range(min(x), (max(x) + 1)):
            coding_space.append(y)
    profile_data = profile_data.reindex(coding_space)

    profile_data.to_csv(
        str(args_dict['coverage']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

"""Get coverage profile for a specific gene"""
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
    args_dict = add_directory(
        args_dict,
        'coverage',
        'individual_plots')

    # Get list of bam files from user input
    files = get_files(
        args_dict['input'],
        [str(args_dict['bam_suffix'])])

    # Get samples user specified
    if args_dict['samples'] != None:
        sample_list = []
        for x in args_dict['samples']:
            for y in files:
                if y in x:
                    sample_list.append(y)
                    break
        files = sample_list

    # Get indices
    print('Generating index for gene...')
    chromosome_index, coordinate_index = get_indices(args_dict, record_type=args_dict['type'], gene_name=args_dict['gene_name'])

    # Perform metagene analysis
    print('Generating coverage profiles for each sample across ' + str(args_dict['gene_name']) + '...')
    func = partial(
        get_coverage,
        chromosome_index = chromosome_index,
        coordinate_index = coordinate_index)
    parallelize(
        func,
        files,
        args_dict,
        mod_workers = True)

    # Compile metrics to plot
    print('Plotting...')
    files = get_files(
        str(args_dict['coverage']) + 'metrics/',
        ['_metrics.txt'])

    file_number = ceil(len(files) / 6)
    file_lists = []

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+16])
        y += 16

    for file_list in file_lists:

        # Plot metrics for each file
        compile_coverage(
            args_dict,
            str(args_dict['coverage']) + 'metrics/',
            file_list,
            chromosome_index,
            coordinate_index,
            'transcript',
            'count',
            'coverage',
            args_dict['experiment'],
            args_dict['coverage'],
            str(args_dict['coverage']) + 'individual_plots/')

    gc.collect()
