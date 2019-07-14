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
    gtf):

    gtf = gtf[(gtf[gtf_type_column] == 'exon') | (gtf[gtf_type_column] == 'CDS')]

    # parse out transcript id to column
    gtf[gtf_transcript_column] = gtf[gtf_annotation_column].str.extract(gtf_transcript_search)

    # only take the essentials
    gtf = gtf[[gtf_transcript_column, gtf_type_column, gtf_sign_column, gtf_leftCoordinate_column, gtf_rightCoordinate_column]]

    # finish formatting
    gtf.columns = ['transcript', 'feature', 'strand', 'left_coordinate', 'right_coordinate']
    gtf = gtf.reset_index(drop=True)

    return gtf

"""Get meta and periodicity indices from GTF"""
def index_gtf(
    args_dict,
    gene_name=None,
    threads=None,
    output=False):

    # Import GTF reference file
    if str(args_dict['gtf']).endswith('.gtf'):
        gtf = pd.read_csv(
            str(args_dict['gtf']),
            sep='\t',
            header=None,
            comment='#',
            low_memory=False)
    else:
        raise Exception('Error: A GTF-formatted file or dataframe was not provided')

    if gene_name != None:
        gene_name = gene_name.replace(' ','')
        output_file = str(gene_name) + '.idx'
        gtf = gtf.loc[gtf[gtf_annotation_column].str.contains(str(gene_name))]
        gtf = gtf.reset_index(drop=True)
    else:
        output_file = 'metagene.idx'

    # Flatten GTF
    if args_dict['gtf'].endswith('LC.gtf') == True:
        gtf_flat = make_index(
            gtf)

    else:
        gtf = edit_gtf(
            gtf,
            longest_transcript=True,
            protein_coding=True,
            truncate_reference=False,
            output=False,
            threads=None)
        gtf_flat = make_index(
            gtf)

    # Get rid of old GTF
    gtf = None
    gc.collect()

    gtf_flat.to_csv(
        str(args_dict['output']) + str(output_file),
        sep = '\t',
        index = False,
        quoting = csv.QUOTE_NONE)

    if output == True:
        return gtf_flat
    else:
        return
