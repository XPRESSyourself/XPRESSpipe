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
from .utils import add_directory
from .parallel import parallelize
from .compile import compile_size_distribution

"""
DESCRIPTION: Generate periodicity maps
"""
def get_peaks(alignments, args_dict):

    #read sam files, may need to convert to bed, get most abundant length and make new file, then perform plastid analysis
    print('')

def make_periodicity(args_dict):

    args_dict = add_directory(args_dict, 'output', 'periodicity')



"""
DESCRIPTION: Generate metagene profiles
"""
def get_profiles(args):

    file, args_dict = args[0], args[1]

    os.system("CollectRnaSeqMetrics -QUIET true -REF_FLAT " + str(args_dict['reference_type']) + " -STRAND_SPECIFICITY NONE -INPUT " + str(args_dict['input']) + str(file) + " -OUTPUT " + str(args_dict['metrics']) + str(file[:-4]) + "_rna_metrics")

def make_metagene(sam_directory, args_dict):

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
