"""XPRESSpipe
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

"""IMPORT DEPENDENCIES
"""
import os
import sys
import pandas as pd
from .utils import add_directory

"""Description:
Purpose of this module is to perform a Ribomap-like isoform quantification for ribosome profiling reads

Prereqs:
Run both HTSeq and Cufflinks (make this an option for quant by providing \"both\")
Quantification should have been performed using a GTF that was truncated and protein-coding-only

Inputs:
counts_table
abundance_table
meta_table -> should include two columns:
    index: sample names as in tables
    1. type

Outputs:
normalized_table

"""


def make_file_dictionary(
        metadata):

    # Get files dictionary
    files = {}

    for x in set(metadata[sample_col].tolist()):

        rpf_name = metadata.loc[(metadata[sample_col] == x) & (metadata[type_col] == 'rpf')].index.tolist()
        rna_name = metadata.loc[(metadata[sample_col] == x) & (metadata[type_col] == 'rna')].index.tolist()

        if len(rpf_name) > 1 or len(rna_name) > 1:
            raise Exception('Only one RPF and one RNA file allowed per sample ID')

        files[x] = {
            'rpf': rpf_name,
            'rna': rna_name}

    return files

def run_riboforms(args):

    sample_id, rpf_file, rna_file, output, args_dict = args[0][0], args[0][1], args[0][2], args[0][3], args[1]

    print('Evaluating the ribosome isoform abundances for ' + str(sample_id))
    # Perform ribosome isoform analysis
    # Loop through each selected BAM

    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'Rriboforms.r'
        + ' ' + str(rpf_file)
        + ' ' + str(rna_file)
        + ' ' + str(output)
        + str(args_dict['log']))

def fetch_files(
        input_bams,
        bam_suffix,
        input_isoforms,
        cufflinks_suffix,
        file):

    # Get RPF file and path
    rpf_file = glob.glob(input_bams + files[key]['rpf'] + '*' + bam_suffix)
    if len(rpf_file) == 1:
        if os.path.isfile(rpf_file[0]):
            rpf_file = rpf_file[0]
        else:
            raise Exception('Failed to find RPF file')
    else:
        raise Exception('Failed to find a single RPF file for sample', key)

    # Get RNA file and path
    rna_file = glob.glob(input_isoforms + files[key]['rpf'] + '*' + cufflinks_suffix)
    if len(rna_file) == 1:
        if os.path.isfile(rna_file[0]):
            rna_file = rna_file[0]
        else:
            raise Exception('Failed to find RNA file')
    else:
        raise Exception('Failed to find a single RNA file for sample', key)

    return rpf_file, rna_file

def __main__(
        args_dict,
        bam_suffix='toTranscriptome.out.bam',
        cufflinks_suffix='_cufflinks_counts'):

    # Set some variables
    sample_info = args_dict['samples']
    input_bams = args_dict['alignments_transcriptome']
    input_isoforms = args_dict['abundances']

    # Step 0: Parse metadata table
    # Requires following fields:
    # 1. sample_id (expected two instances each)
    # 2. Sequencing type (RPF or RNA)
    metadata = pd.read_csv(
        sample_info,
        sep='\t',
        index_col=0)

    if len(metadata.columns.tolist) > 2:
        raise Exception('Expected two columns only in metadata table')

    # Get seq type and sample id column names
    type_col = ''
    for x in metadata.columns.tolist():
        if 'rpf' in [x.lower() for x in metadata[x].tolist()] \
        and 'rna' in [x.lower() for x in metadata[x].tolist()]:
            type_col = x
        else:
            sample_col = x

    metadata[type_col] = [x.lower() for x in metadata[type_col].tolist()]
    files = make_file_dictionary(
        metadata=metadata)

    # Add output directories
    args_dict = add_directory(
        args_dict,
        'output',
        'ribo_abundances')
    output = args_dict['ribo_abundances']

    # Add lists of dictionaries to process
    file_lists = []

    for key in files.keys():

        rpf_file, rna_file = fetch_files(
            input_bams,
            bam_suffix,
            input_isoforms,
            cufflinks_suffix,
            file)
        file_lists.append([key, rpf_file, rna_file, output])

    # Perform metagene analysis
    parallelize(
        run_riboforms,
        file_lists,
        args_dict,
        mod_workers = True)

    # Combine abundances from RPFs and RNAs to master table for output
