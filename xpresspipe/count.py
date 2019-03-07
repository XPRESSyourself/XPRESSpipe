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
from xpresstools import count_table, batch_normalize
from .utils import get_files, add_directory
from .parallel import parallelize

"""
DESCRIPTION: Convert a sorted sam file to a bam file
"""
def sam2bam(path, file):

    os.system('samtools view -h -S -b ' + str(path) + str(file) + ' > ' + str(path) + str(file[:-4]) + '.bam')
    os.system("samtools index " + str(path) + str(file[:-4]) + '.bam')

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
    print('bedtools bamtobed -i ' + str(args_dict['input']) + str(file[:-4]) + '.bam > ' + str(args_dict['bed_files']) + str(file[:-4]) + '.bed')
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
    os.system('bamCoverage -b ' + str(args_dict['input']) + str(file[:-4]) + '.bam -o ' + str(args_dict['bigwig_files']) + str(file[:-4]) + '.bw')

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

    #Determine reference file type
    if args_dict['count_coding'] == True and args_dict['truncate'] == True:
        transcript_type = 'transcripts_coding_truncated'

    elif args_dict['count_coding'] == True and args_dict['truncate'] == False:
        transcript_type = 'transcripts_coding'

    elif args_dict['count_coding'] == False and args_dict['truncate'] == True:
        raise Exception('A truncated transcript file that is not coding-only transcripts is not currently supported')

    elif args_dict['count_coding'] == False and args_dict['truncate'] == False:
        transcript_type = 'transcripts'

    else:
        raise Exception('Something went wrong in determining transcript reference file type')

    #Count
    if args_dict['type'] == 'PE':
        os.system('htseq-count -r pos -s no ' + str(args_dict['input']) + str(file) + ' ' + str(args_dict['reference']) + str(transcript_type)+ '.gtf > ' + str(args_dict['counts']) + str(file[:-4]) + '.tsv')
    else:
        os.system('htseq-count -s no ' + str(args_dict['input']) + str(file) + ' ' + str(args_dict['reference']) + str(transcript_type)+ '.gtf > ' + str(args_dict['counts']) + str(file[:-4]) + '.tsv')

def count_reads(args_dict):

    #Add output directories
    args_dict = add_directory(args_dict, 'output', 'counts')

    #Get list of files to count based on acceptable file types
    files = get_files(args_dict['input'], ['.sam'])

    #Count aligned RNAseq reads
    parallize(count_file, files, args_dict)

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
    c = 0
    for x in files:
        files[c] = str(args_dict['input']) + str(x)
        c += 1

    #Create and output collated count table
    count_table = count_table(files)
    if 'experiment' in args_dict:
        count_table.to_csv(str(args_dict['experiment']) + 'counts_table.csv')
    else:
        cdt = datetime.datetime.now()
        count_table.to_csv(str(cdt.year) + '_' + str(cdt.month) + '_' + str(cdt.day) + '_' + str(cdt.hour) + 'h_' + str(cdt.minute) + 'm_' + str(cdt.second) + 's_counts_table.csv')

"""
DESCRIPTION: Run normalization of count dataframe
"""
def run_normalization(args_dict):
    if 'method' in args_dict:
        #RPM normalization
        if args_dict['method'].upper() == 'RPM':
            type = 'rpm'
            df = pd.read_csv(str(args_dict['data']), sep=',', header=0, index_col=0, comment='#', low_memory=False)
            df = rpm(df)
            df.to_csv(str(args_dict['data'][:-4]) + '_' + str(type) + 'Normalized.csv', sep=',')
        #RPKM or FPKM normalization
        elif args_dict['method'].upper() == 'RPKM' or args_dict['type'].upper() == 'FPKM':
            type = 'r_fpkm'
            df = pd.read_csv(str(args_dict['data']), sep=',', header=0, index_col=0, comment='#', low_memory=False)
            df = r_fpkm(df)
            df.to_csv(str(args_dict['data'][:-4]) + '_' + str(type) + 'Normalized.csv', sep=',')
        #Log normalization
        elif args_dict['method'].upper() == 'LOG':
            type = 'log'
            df = pd.read_csv(str(args_dict['data']), sep=',', header=0, index_col=0, comment='#', low_memory=False)
            df = log_scale(df, log_base=10)
            df.to_csv(str(args_dict['data'][:-4]) + '_' + str(type) + 'Normalized.csv', sep=',')
        else:
            raise Exception('Unknown \"method\" argument provided')
    #Run in batch normalization
        if 'batch' in args_dict:
            batch_normalize(str(args_dict['data'][:-4]) + '_' + str(type) + 'Normalized.csv', str(args_dict['batch']), str(args_dict['output']), input_sep=',', batch_sep=',')
    else:
        if 'batch' in args_dict:
            batch_normalize(str(args_dict['data']), str(args_dict['batch']), str(args_dict['output']), input_sep=',', batch_sep='\t')
