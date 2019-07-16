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
import gc
import csv
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES"""
from .gtfModify import edit_gtf

"""GLOBALS"""
gtf_chr_column = 0
gtf_type_column = 2
gtf_leftCoordinate_column = 3
gtf_rightCoordinate_column = 4
gtf_sign_column = 6
gtf_annotation_column = 8
gtf_transcript_search = r'transcript_id \"(.*?)\"; '
gtf_transcript_column = 9

"""Flatten GTF dataframe"""
def make_index(
    args_dict,
    gtf_file,
    output):

    os.system(
        'Rscript'
        + ' ' + str(args_dict['path']) + 'RbuildIndex.r'
        + ' ' + str(gtf_file)
        + ' ' + str(output)
        + str(args_dict['log']))

"""Make feature index for a gene"""
def make_features(
    gtf,
    output):

    # Pull out relevant records
    gtf = gtf[(gtf[gtf_type_column] == 'exon') | (gtf[gtf_type_column] == 'CDS') | (gtf[gtf_type_column] == 'five_prime_utr') | (gtf[gtf_type_column] == 'three_prime_utr')]

    # Parse out transcript id to column
    gtf[gtf_transcript_column] = gtf[gtf_annotation_column].str.extract(gtf_transcript_search)
    gtf['length'] = gtf[gtf_rightCoordinate_column] - gtf[gtf_leftCoordinate_column] + 1

    # Only take the essentials
    gtf = gtf[[gtf_transcript_column, gtf_type_column, 'length']]

    # Finish formatting
    gtf.columns = ['transcript', 'feature', 'length']
    gtf = gtf.reset_index(drop=True)

    # Output for use by modules
    gtf.to_csv(
        str(output),
        sep = '\t',
        index = False,
        quoting = csv.QUOTE_NONE)

"""Get meta and periodicity indices from GTF"""
def index_gtf(
    args_dict,
    gene_name=None):

    # Check for GTF file input
    if str(args_dict['gtf']).endswith('.gtf'):
        print('Generating index for genes...')
    else:
        raise Exception('Error: A GTF-formatted file or dataframe was not provided')

    # Gene coverage index creation
    if gene_name != None:

        # Get file names and clean up inputs
        gene_name = gene_name.replace(' ','')
        gene_gtf = str(gene_name) + '.gtf'
        output_file = str(gene_name) + '.idx'

        # Import GTF and get only records for gene of interest
        gtf = pd.read_csv(
            str(args_dict['gtf']),
            sep='\t',
            header=None,
            comment='#',
            low_memory=False)

        gtf = gtf.loc[gtf[gtf_annotation_column].str.contains('\"' + str(gene_name) + '\"')]
        gtf = gtf.reset_index(drop=True)

        # Get canonical, longest only GTF for the given gene
        gtf = edit_gtf(
            gtf,
            longest_transcript=True,
            protein_coding=False,
            truncate_reference=False,
            output=False,
            threads=None)

        # Save GTF file import by other functions
        gtf.to_csv(
            str(args_dict['output']) + str(gene_gtf),
            sep = '\t',
            header = None,
            index = False,
            quoting = csv.QUOTE_NONE)

        # Make features index for geneCoverage plotting of exons and CDS start/stop
        make_features(
            gtf,
            str(args_dict['output']) + str(gene_name) + '.fts')
        gtf = None

        # Make index for gene to calculate coverage
        make_index(
            args_dict,
            str(args_dict['output']) + str(gene_gtf),
            str(args_dict['output']) + str(gene_name) + '.idx')

        # Remove intermediate GTF file
        os.system(
            'rm'
            + ' ' + str(args_dict['output']) + str(gene_gtf))

    else:
        # Make index of all transcripts
        make_index(
            args_dict,
            args_dict['gtf'],
            str(args_dict['output']) + 'transcripts.idx')
