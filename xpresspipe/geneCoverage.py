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
from .processBAM import read_bam, bam_sample
from .quality import get_indices

plots_per_page = 8
chromosome_position = 2
leftCoordinate_position = 3
read_position = 9

"""
func: Find right-most coordinate for each read
@param bam: BAM-format pandas dataframe
@return: multi-dimensional array of chromosome, start, and stop coordinates
"""
def cov_coordinates(
    bam):

    # Map middle point of each read as left-most position plus half of read length
    bam['right_coordinate'] = bam[leftCoordinate_position] + (bam[read_position].str.len()).astype('int64')

    # Return as array
    mid_coordinates = bam[[2,3,'right_coordinate']]
    return mid_coordinates.values.tolist()

"""
func: Find if a given coordinate falls in a gene exon range
@param position: genome coordinate
@param coordinates: a list of lists of start and stop coordinates for each exon of a gene
@return: Position if falls in exon coordinate range, or None if does not map
"""
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

"""Get meta profile for bam file
- Assumes input bam index only grabbed from the chromosome of interest
- coordinate_index has already been selected based on chromosome
- Counts the coordinate for every nucleotide of a read that falls within an exon for the given gene
@param aligned_reads_index: Multidimensional array with chromosome number (positional), leftmost position, and rightmost position for each record in BAM
@param coordinate_index: Gene index array
@return: Metagene profile dataframe with mapped counts
"""
def get_cov_profile(
    aligned_reads_index,
    coordinate_index):

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
- Parse through BAM file for plausible reads falling within specified gene region
- Generate a coverage profile dataframe based on this information
@param args: args iterator for each file in incoming list with file name and args_dict
@param chromosome_index: chromosome index info for gene of interest
@param coordinate_index: coordindate index info for gene of interest
@return: None, saves file as metrics table
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

    # Check that there is only one gene record for gene coverage
    # Get plausible gene locations in BAM file for gene record of interest
    if len(coordinate_index) == 1 and len(coordinate_index[0]) == 1 and len(chr) == 1:
        bam = bam.loc[(bam[chromosome_position] == chr[0]) & (bam[leftCoordinate_position] >= (coordinate_index[0][0][0] - coordinate_index[0][0][4])) & (bam[leftCoordinate_position] <= (coordinate_index[0][0][1] + coordinate_index[0][0][4]))]
    else:
        raise Exception('A single gene record was not selected.')

    # Get meta-data for bam file relevant to gene region
    coordinates = cov_coordinates(bam)

    bam = None # Some clean-up
    del bam

    # Get profile
    profile_data = get_cov_profile(
        coordinates,
        coordinate_index)

    coding_space = []
    for x in coordinate_index[0][0][3]:
        for y in range(min(x), (max(x) + 1)):
            coding_space.append(y)
    profile_data = profile_data.reindex(coding_space)

    profile_data.to_csv(
        str(args_dict['coverage']) + 'metrics/' + str(file)[:-4] + '_metrics.txt',
        sep='\t')

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
    chromosome_index, coordinate_index = get_indices(
        args_dict,
        record_type=args_dict['type'],
        gene_name=args_dict['gene_name'],
        threads=1)

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
            args_dict['type'],
            args_dict['sample_names'],
            coordinate_index[0][0][2],
            'coverage' + str(z),
            args_dict['experiment'],
            args_dict['coverage'],
            args_dict['plot_color'])

        z += 1

    gc.collect()
