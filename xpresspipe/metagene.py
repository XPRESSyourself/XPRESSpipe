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
from .buildIndex import index_gtf

transcript_length = 'l_tr'
utr5_length = 'l_utr5'
utr3_length = 'l_utr3'
cds_length = 'l_cds'

def roundup(value):
    return int(ceil(value))

def finish_metagene(
        args,
        remove_outliers=True):

    file, args_dict = args[0], args[1]

    # Read in parsed bam
    bam = pd.read_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '.metaposit',
        sep = '\t')

    # Read in transcript exon space dictionary
    dict = pd.read_csv(
        str(args_dict['output']) + 'transcripts.idx',
        sep = '\t',
        index_col=0)

    # Make transcript / length dictionary
    reference = pd.Series(dict[transcript_length].values, index=dict['transcript']).to_dict()

    # Map exon space to each bam record based on its transcript ID
    bam['total_length']= bam['seqnames'].map(reference)

    # Get UTR offset amounts and map to transcript IDs
    if str(args_dict['feature_type']).lower() == 'cds':
        utr5_correct = pd.Series(dict[utr5_length].values, index=dict['transcript']).to_dict()
        utr3_correct = pd.Series(dict[utr3_length].values, index=dict['transcript']).to_dict()
        cds_correct = pd.Series(dict[cds_length].values, index=dict['transcript']).to_dict()

        bam['utr5_offset'] = bam['seqnames'].map(utr5_correct)
        bam['utr3_offset'] = bam['seqnames'].map(utr3_correct)
        bam['cds_offset'] = bam['seqnames'].map(cds_correct)

    dict = None
    bam = bam.dropna()

    # Calculate meta-position
    if str(args_dict['feature_type']).lower() == 'cds':
        # Remove records where meta_position is before CDS start
        bam = bam.loc[(bam['meta_position'] - bam['utr5_offset']) >= 0]

        # Remove records where meta_position is after CDS stop
        bam = bam.loc[(bam['meta_position'] - bam['utr5_offset'] - bam['cds_offset']) < 0]

        # Compensate for UTRs in meta-calculations across CDS
        # Those reads mapping outside of the CDS region will be < 1 or > 100
        # When data is collated, will only take values between 1 and 100
        bam['meta_distance'] = (bam['meta_position'] - bam['utr5_offset']) / (bam['cds_offset']) * 100
    else:
        bam['meta_distance'] = (bam['meta_position']) / (bam['total_length']) * 100

    # Set a nice number with no remainder
    bam['meta_distance'] = bam['meta_distance'].apply(roundup)

    # Get rid of outlier-expressors
    if remove_outliers == True:
        # Will count the number of members in each group, will be the new meta_distance column of this dataframe
        bam_outliers = bam.groupby('seqnames').count().sort_values('meta_distance')
        del bam_outliers.index.name # This is seqnames label

        removal_amount = int(ceil(bam_outliers.shape[0] * 0.005))
        bam_outliers = bam_outliers.iloc[removal_amount:-removal_amount]

        # Take list of transcripts that made the cut and only grab those out of the starting BAM dataframe
        bam = bam[bam['seqnames'].isin(bam_outliers.index.tolist())]

    # Compile meta position statistics
    profile = pd.DataFrame()
    profile['metacount'] = [0 for x in range(1,101)]
    profile.index = range(1,101)

    for x in range(1,101):
        # Automatically filters out reads mapped outside of CDS as they will be < 0 or > 1
        # Go through each position in the meta-transcript
        # Get number of members for each transcript at that position
        # Take the mean of those values -- reduces noise of one transcript with 100x more reads than the next highest
        # Outlier expressors were previously parsed out
        profile.loc[x] = bam[bam['meta_distance'] == x].groupby('seqnames').count().mean()[0]

    # Set mean to sum of meta-positions = 100
    profile['metacount'] = profile['metacount'] / profile['metacount'].mean()

    # Export metrics
    profile['representative transcript'] = profile.index
    profile.to_csv(
        str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '_metrics.txt',
        sep='\t')

    # Clean up
    bam = None
    reference = None
    dict = None
    profile = None
    gc.collect()
    os.system(
        'rm'
        + ' ' + str(args_dict['metagene']) + 'metrics/' + str(file).split('.')[0] + '.metaposit')

def run_metagene(args):

    file, args_dict = args[0], args[1]
    file = '\"' + str(file) + '\"'

    print('Evaluating the metagene profile for ' + str(file))
    # Perform metagene analysis
    # Loop through each selected BAM
    # args[1] = Path to BAM files
    # args[2] = List of BAM files
    # args[3] = Index file with path
    # args[4] = Output file path
    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'Rmetagene.r'
        + ' ' + str(args_dict['input'])
        + ' ' + str(file)
        + ' ' + str(args_dict['metagene']) + 'metrics/'
        + str(args_dict['log']))

"""Manager for running metagene summary plotting"""
def make_metagene(
        args_dict,
        files):

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

    # Perform metagene analysis
    print('\nGetting positions of reads...')
    parallelize(
        run_metagene,
        files,
        args_dict,
        mod_workers = True)

    print('\nConverting positions to metapositions...')
    parallelize(
        finish_metagene,
        files,
        args_dict,
        mod_workers = True)

    # Compile metrics to plot
    print('\nPlotting...')
    files = get_files(
        str(args_dict['metagene']) + 'metrics/',
        ['_metrics.txt'])

    file_number = ceil(len(files) / 6)
    file_lists = []

    y = 0
    for x in range(file_number):
        file_lists.append(files[y:y+6])
        y += 6

    z = 1
    for file_list in file_lists:

        # Plot metrics for each file
        compile_matrix_metrics(
            str(args_dict['metagene']) + 'metrics/',
            file_list,
            'representative transcript',
            'metacount',
            str(args_dict['feature_type']) + '_metagene_' + str(z),
            args_dict['experiment'],
            args_dict['metagene'])
        z += 1

    gc.collect()
