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
import re
import csv
import gc
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES
"""
from .utils import get_files
from .parallel import parallelize

gtf_type_column = 2
gtf_annotation_column = 8
search_type = 'transcript'
parse_type = 'rrna'

def generate_bed(
        gtf_file):
    """Generate a BED file of rRNA sequences for depletion from genome-aligned
    BAM file.
    """

    # Set up BED generation command
    gtf = pd.read_csv(
        str(gtf_file),
        sep = '\t',
        header = None,
        comment = '#',
        low_memory = False)

    # Remove records that map to rRNA
    gtf_rrna = gtf[gtf[gtf_annotation_column].str.contains(
        parse_type, flags=re.IGNORECASE
        )]
    gtf_rrna = gtf_rrna.loc[gtf_rrna[gtf_type_column] == search_type]
    gtf_rrna = gtf_rrna[[0,3,4]]

    bed_file = gtf_file[:-4] + '_rrna.bed'

    gtf_rrna.to_csv(
        str(bed_file),
        sep = '\t',
        header = None,
        index = False,
        quoting = csv.QUOTE_NONE)

    gtf = None # Garbage management
    gtf_rrna = None
    gc.collect()

    return bed_file

def run_intersect(
        args):
    """Use rRNA BED file to extract the inverse intersection of alignments
    args[0] = alignment file directory
    args[1] = bam_file of interest
    args[2] = bed_file for reference
    """

    file_iterator, args_dict = args[0], args[1]
    dir, bam_file, bed_file = file_iterator[0], file_iterator[1], file_iterator[2]

    # Set up intersection command
    cmd = (
        'bedtools intersect -abam ' + dir + bam_file
        + ' -b ' + bed_file
        + ' -v > ' + dir + bam_file + '.depl.bam'
    )

    # Run
    os.system(cmd)

    os.system(
        'mv '
        + dir + bam_file + '.depl.bam'
        + ' ' + dir + bam_file
    )

def genomic_depletion(
        args_dict,
        dir='alignments_coordinates',
        suffix='_Aligned.sort.bam'):
    """Remove rRNA records from genome-aligned BAM file

    Usage:
    - Access gtf file via:
        str(args_dict['reference']) + 'transcripts.gtf'
    - Access genome-aligned files via:
        str(args_dict['alignments_coordinates'])
    * file_list should already have path appended
    """

    gtf_file = str(args_dict['reference']) + 'transcripts.gtf'
    bed_file = generate_bed(
        gtf_file=gtf_file)

    # Get list of files to count based on acceptable file types
    file_list = get_files(
        args_dict[dir],
        [str(suffix)])

    file_iterator = []
    for file in file_list:

        file_iterator.append([args_dict[dir], file, bed_file])

    parallelize(
        run_intersect,
        file_iterator,
        args_dict,
        mod_workers = 'all')
