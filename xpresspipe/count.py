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

"""
IMPORT DEPENDENCIES
"""
import os, sys
import datetime
import pandas as pd
from xpresstools import count_table, batch_normalize, rpm, r_fpkm, log_scale
from .utils import get_files, add_directory
from .parallel import parallelize

"""
DESCRIPTION: Convert a sorted sam file to a bam file
"""
def sam2bam(path, file):

    os.system('samtools view -h -S -b ' + str(path) + str(file) + ' > ' + str(path) + str(file[:-4]) + '.bam')
    os.system("samtools index " + str(path) + str(file[:-4]) + '.bam' + str(args_dict['log']))

"""
DESCRIPTION: Convert sorted sam files in directory to bed files
"""
def bed_convert(args):

    file, args_dict = args[0], args[1]

    #Ensure input file is properly formatted as a sorted and indexed BAM file
    if file.endswith('.sam'):
        sam2bam(args_dict['input'], file)
    elif file.endswith('.bam'):
        pass
    else:
        raise Exception('Incorrect input file')

    #Convert BAM to BED
    os.system('bedtools bamtobed -i ' + str(args_dict['input']) + str(file[:-4]) + '.bam > ' + str(args_dict['bed_files']) + str(file[:-4]) + '.bed')

def create_bed(args_dict):

    #Add output directories
    args_dict = add_directory(args_dict, 'output', 'bed_files')

    #Get list of files to convert based on acceptable file types
    files = get_files(args_dict['input'], ['.sam', '.bam'])

    #Convert aligned RNAseq reads to BED files
    parallelize(bed_convert, files, args_dict)

    return args_dict

"""
DESCRIPTION: Convert sorted sam files in directory to bigwig files
"""
def bigwig_convert(args):

    file, args_dict = args[0], args[1]

    #Ensure input file is properly formatted as a sorted and indexed BAM file
    if file.endswith('.sam'):
        sam2bam(args_dict['input'], file)
    elif file.endswith('.bam'):
        pass
    else:
        raise Exception('Incorrect input file')

    #Convert BAM to bigwig
    os.system('bamCoverage -b ' + str(args_dict['input']) + str(file[:-4]) + '.bam -o ' + str(args_dict['bigwig_files']) + str(file[:-4]) + '.bw' + str(args_dict['log']))

def create_bigwig(args_dict):

    #Add output directories
    args_dict = add_directory(args_dict, 'output', 'bigwig_files')

    #Get list of files to convert based on acceptable file types
    files = get_files(args_dict['input'], ['.sam', '.bam'])

    #Convert aligned RNAseq reads to bigwig files
    parallelize(bigwig_convert, files, args_dict)

    return args_dict

"""
DESCRIPTION: Compile counts tables from HTseq output files
"""
def count_file(args):

    file, args_dict = args[0], args[1]

    #Count
    os.system('htseq-count -q -m intersection-nonempty -t exon -r pos -s no ' + str(args_dict['input']) + str(file) + ' ' + str(args_dict['gtf_type']) + ' > ' + str(args_dict['counts']) + str(file[:-4]) + '.tsv')

def count_reads(args_dict):

    #Add output directories
    args_dict = add_directory(args_dict, 'output', 'counts')

    #Get list of files to count based on acceptable file types
    files = get_files(args_dict['input'], ['.sam'])

    #Count aligned RNAseq reads
    parallelize(count_file, files, args_dict, mod_workers=True)

    return args_dict

"""
DESCRIPTION: Take directory of single counts files and collate into single table
"""
def collect_counts(args_dict):

    #Add output directories
    args_dict = add_directory(args_dict, 'output', 'counts')

    #Get list of files to count based on acceptable file types
    files = get_files(args_dict['input'], ['.tsv'])

    #Append path to file list
    count_files = []
    for x in files:
        count_files.append(str(args_dict['input']) + str(x))

    #Create and output collated count table
    counts = count_table(count_files)

    if 'experiment' in args_dict and args_dict['experiment'] != None:
        counts.to_csv(str(args_dict['counts']) + str(args_dict['experiment']) + '_counts_table.tsv', sep='\t')
    else:
        cdt = datetime.datetime.now()
        counts.to_csv(str(args_dict['counts']) + str(cdt.year) + '_' + str(cdt.month) + '_' + str(cdt.day) + '_' + str(cdt.hour) + 'h_' + str(cdt.minute) + 'm_' + str(cdt.second) + 's_counts_table.tsv', sep='\t')

"""
DESCRIPTION: Run normalization of count dataframe
"""
def run_normalization(args_dict, sep='\t'):

    #Run sample normalization
    if 'method' in args_dict and args_dict['method'] != None:
        #RPM normalization
        if args_dict['method'].upper() == 'RPM':
            type = 'rpm'
            df = pd.read_csv(str(args_dict['input']), sep=sep, header=0, index_col=0, comment='#', low_memory=False)
            df = rpm(df)
            df.to_csv(str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv', sep='\t')
        #RPKM or FPKM normalization
        elif args_dict['method'].upper() == 'RPKM' or args_dict['method'].upper() == 'FPKM':
            if args_dict['gtf'] == None:
                raise Exception('A GTF reference file is required for RPKM and FPKM normalization')
            type = 'r_fpkm'
            df = pd.read_csv(str(args_dict['input']), sep=sep, header=0, index_col=0, comment='#', low_memory=False)
            df = r_fpkm(df, args_dict['gtf'])
            df.to_csv(str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv', sep='\t')
        #Log normalization
        elif args_dict['method'].upper() == 'LOG':
            type = 'log'
            df = pd.read_csv(str(args_dict['input']), sep=sep, header=0, index_col=0, comment='#', low_memory=False)
            df = log_scale(df, log_base=10)
            df.to_csv(str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv', sep='\t')
        else:
            raise Exception('Unknown \"method\" argument provided')

    #Run in batch normalization
        if 'batch' in args_dict and args_dict['batch'] != None:
            batch_normalize(str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.csv', str(args_dict['batch']))
    else:
        if 'batch' in args_dict and args_dict['batch'] != None:
            batch_normalize(str(args_dict['input']), str(args_dict['batch']))
