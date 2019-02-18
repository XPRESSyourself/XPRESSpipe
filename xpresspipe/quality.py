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
import pandas as pd
from statistics import mode
from .utils import add_directory
from .parallel import parallelize
from .compile import compile_size_distribution

"""
DESCRIPTION: Generate periodicity maps
"""
def get_peaks(args):

    file, args_dict = args[0], args[1]

    #read sorted, unique only sam files, get most abundant length and make new file
    df = pd.read_csv(str(args_dict['input']) + str(file), sep='\t', comment='@' header=None, usecols = list(range(0, 15)))
    df[15]  = df[9].str.len()
    df_mode = df.loc[df[15] == mode(list(df[15]))]
    df_mode.to_csv(str(args_dict['metrics']) + str(file[:-4]) + 'mode.sam')

    #Perform plastid analysis
    os.system('metagene count -q ' + str(args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]) + 'metagene_reference_' + str(args_dict['downstream']) + '_rois.txt ' + str(args_dict['metrics']) + str(file[:-4]) + '_periodicity --count_files ' + str(args_dict['metrics']) + str(file) + ' --fiveprime --offset 14 --normalize_over 30 200 --min_counts 50 --cmap Blues --title ' + file[:-4])

def make_periodicity(args_dict):

    args_dict = add_directory(args_dict, 'output', 'periodicity')
    args_dict = add_directory(args_dict, 'periodicity', 'metrics')

    #Get list of sam files from user input
    files = get_files(args_dict['input'], ['.sam'])

    #Add metagene reference to folder
    os.system('metagene generate ' + str(args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]) + 'metagene_reference --landmark ' + str(args_dict['landmark']) + ' --annotation_files ' + str(args_dict['gtf']) + ' --downstream ' + str(args_dict['downstream']))

    #Perform metagene analysis
    parallize(get_peaks, files, args_dict)

    #Compile images
    files = get_files(args_dict['metrics'], ['_periodicity.txt'])
    compile_size_distribution(args_dict, args_dict['metrics'], files, 'metagene_average', str(int(args_dict['downstream'] - 1)), 'meta_position', 'normalized_coverage (all_reads)', 'periodicity', args_dict['experiment'])

    #Clean mode sams
    os.system('rm ' + str(args_dict['metrics']) + '*_mode.sam')

"""
DESCRIPTION: Generate metagene profiles
"""
def get_profiles(args):

    file, args_dict = args[0], args[1]

    os.system("CollectRnaSeqMetrics -QUIET true -REF_FLAT " + str(args_dict['reference_type']) + " -STRAND_SPECIFICITY NONE -INPUT " + str(args_dict['input']) + str(file) + " -OUTPUT " + str(args_dict['metrics']) + str(file[:-4]) + "_rna_metrics")

def make_metagene(args_dict):

    #Add output directory to output for metagene profiles
    args_dict = add_directory(args_dict, 'output', 'metagene')
    args_dict = add_directory(args_dict, 'metagene', 'metrics')

    #Get refFlat file type
    if args_dict['reference_type'].upper() == 'DEFAULT':
        args_dict['reference_type'] = str(args_dict['reference']) + 'transcripts_refFlat.txt'
    elif args_dict['reference_type'].upper() == 'CODING':
        args_dict['reference_type'] = str(args_dict['reference']) + 'transcripts_coding_refFlat.txt'
    elif args_dict['reference_type'].upper() == 'CODING_TRUNCATED':
        args_dict['reference_type'] = str(args_dict['reference']) + 'transcripts_coding_truncated_refFlat.txt'
    else:
        raise Exception('Invalid reference flat file type provided')

    #Get list of sam files from user input
    files = get_files(args_dict['input'], ['.sam'])

    #Perform metagene analysis
    parallize(get_profiles, files, args_dict)

    #Compile images
    files = get_files(args_dict['metrics'], ['_rna_metrics.txt'])
    compile_size_distribution(args_dict, args_dict['metrics'], files, 'normalized_position', None, 'meta_position', 'normalized_coverage (all_reads)', 'metagene', args_dict['experiment'])

"""
DESCRIPTION: Generate read distribution profiles
"""
def run_fastqc(args):

    file, args_dict = args[0], args[1]

    os.system("fastqc -q " + str(args_dict['input']) + str(file) + " -o " + str(args_dict['fastqc_output']))

def make_readDistributions(args_dict):

    args_dict = add_directory(args_dict, 'output', 'read_distributions')
    args_dict = add_directory(args_dict, 'read_distributions', 'fastqc_output')

    #Get FASTQC file list and unzip
    files = get_files(args_dict['input'], ['.fastq', '.fq', '.txt'])

    #Perform fastqc on each file and unzip output
    parallize(run_fastqc, files, args_dict)

    files = get_files(args_dict['fastqc_output'], ['.zip'])
    for file in files:
        if file.endswith('.zip'):
            os.system('unzip -q ' + str(args_dict['fastqc_output']) + str(file) + ' -d ' + str(args_dict['fastqc_output']))

    #Compile read distributions
    files = get_files(args_dict['fastqc_output'], ['.zip'])
    c = 0
    for x in files:
        files[c] = x[:-4] + '/fastqc_data.txt'
        c += 1
    compile_size_distribution(args_dict, args_dict['fastqc_output'], files, '#Length', '>>END_MODULE', 'position', 'read size (bp)', 'fastqc', args_dict['experiment'])
