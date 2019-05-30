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

"""IMPORT DEPENDENCIES"""
import os
import sys
import datetime
import pandas as pd
from xpressplot import count_table

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, get_directories, add_directory
from .parallel import parallelize

__path__ = str(os.path.dirname(os.path.realpath(__file__))) + '/'

"""Parse cufflinks table for FPKM info"""
def parse_table(
    directory,
    output,
    file):

    # file will be at str(directory) + str(file)
    # output will be to str(output) as directory_name.tsv from second to last /
    output_file = str(output) + str(directory)[:-1].split('/')[-1] + '_fpkm_counts.tsv'

    # Read in cufflinks quantification table
    cufflinks_table = pd.read_csv(
        str(directory) + str(file),
        sep = '\t')

    # Get relevant information
    cufflinks_out = cufflinks_table[['tracking_id','FPKM']]
    cufflinks_out.columns = [0,1]

    # Output table to specified directory with XPRESSpipe parsable file name
    cufflinks_out.to_csv(
        str(output_file),
        sep = '\t',
        index = False,
        header = False)

    return 0

"""Create counts tables from HTseq"""
def count_file_htseq(
    args):

    file, args_dict = args[0], args[1] # Parse args

    # Count
    os.system(
        'htseq-count'
        + ' -q'
        + ' -f bam'
        + ' -m intersection-nonempty'
        + ' -t exon'
        + ' -i gene_id'
        + ' -r pos'
        + ' -s no'
        + ' ' + str(args_dict['input']) + str(file)
        + ' ' + str(args_dict['gtf'])
        + ' > ' + str(args_dict['counts']) + str(file[:-4]) + '.tsv')

"""Create counts tables from cufflinks"""
def count_file_cufflinks(
    args):

    file, args_dict = args[0], args[1] # Parse args

    # Add output directories
    args_dict = add_directory(
        args_dict,
        'counts',
        str(file[:-4]) + '_cufflinks_counts')

    # Count
    os.system(
        str(__path__) + 'cufflinks'
        + ' ' + str(args_dict['input']) + str(file)
        + ' --output-dir ' + str(args_dict['counts']) + str(file[:-4]) + '_cufflinks_counts'
        + ' --GTF ' + str(args_dict['reference']) + 'transcripts.gtf'
        + ' --num-threads 1'
        + str(args_dict['log']))

"""Run count reads manager"""
def count_reads(
    args_dict,
    suffix='_Aligned.sort.bam'):

    # Add output directories
    args_dict = add_directory(
        args_dict,
        'output',
        'counts')

    if 'deduplicate' in args_dict and args_dict['deduplicate'] == True:
        suffix = '_dedupRemoved.bam'

    if 'bam_suffix' in args_dict:
        suffix = args_dict['bam_suffix']

    # Get list of files to count based on acceptable file types
    files = get_files(
        args_dict['input'],
        [str(suffix)])

    # Count aligned RNAseq reads
    if args_dict['quantification_method'] == 'cufflinks':
        parallelize(
        count_file_cufflinks,
        files,
        args_dict,
        mod_workers = True)
    else:
        parallelize(
        count_file_htseq,
        files,
        args_dict,
        mod_workers = True)

    return args_dict

"""Take directory of single counts files and collate into single table"""
def collect_counts(
    args_dict):

    # Add output directories
    if 'counts' not in args_dict:
        args_dict = add_directory(
            args_dict,
            'output',
            'counts')

    # Get list of files to count based on acceptable file types
    if args_dict['quantification_method'] == 'cufflinks':

        directories = get_directories(
            args_dict['input'],
            ['_cufflinks_counts'])

        for dir in directories:
            parse_table(
                dir,
                args_dict['input'],
                'genes.fpkm_tracking')

        files = get_files(
            args_dict['input'],
            ['fpkm_counts.tsv'])

    else:
        files = get_files(
            args_dict['input'],
            ['.tsv'])

    # Append path to file list
    count_files = []
    for x in files:
        count_files.append(str(args_dict['input']) + str(x))

    # Create and output collated count table
    counts = count_table(count_files)

    # Output data
    if args_dict['quantification_method'] == 'cufflinks':
        count_suffix = '_cufflinksFPKM_table.tsv'
    else:
        count_suffix = '_count_table.tsv'

    if 'experiment' in args_dict and args_dict['experiment'] != None:
        counts.to_csv(
            str(args_dict['counts']) + str(args_dict['experiment']) + str(count_suffix),
            sep = '\t')
    else:
        cdt = datetime.datetime.now()
        counts.to_csv(
            str(args_dict['counts']) + str(cdt.year) + '_' + str(cdt.month) + '_' + str(cdt.day) \
            + '_' + str(cdt.hour) + 'h_' + str(cdt.minute) + 'm_' + str(cdt.second) + 's' + str(count_suffix),
            sep = '\t')
