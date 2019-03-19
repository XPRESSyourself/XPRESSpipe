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
import csv
import pandas as pd
from statistics import mode
from .parallel import parallelize
from .compile import compile_size_distribution
from .utils import add_directory, get_files

"""
DESCRIPTION: Generate periodicity maps
"""
def get_peaks(args):

    file, args_dict = args[0], args[1]

    #Get SAM headers
    with open(str(args_dict['input']) + str(file)) as comment:
        reader = csv.reader(comment)
        comment1 = next(reader)
        comment2 = next(reader)
        comment3 = next(reader)
        comment4 = next(reader)

    #read sorted, unique only sam files, get most abundant length and make new file
    df = pd.read_csv(str(args_dict['input']) + str(file), sep='\t', skiprows=4, header=None, usecols = list(range(0, 16)))
    df[16]  = df[9].str.len()
    df_mode = df.loc[df[16] == mode(list(df[16]))]
    df_mode = df_mode.drop(labels=16, axis=1)

    #Save SAM file with SAM header
    f = open(str(args_dict['metrics']) + str(file[:-4]) + '_mode.txt', 'a')
    f.write(str(''.join(comment1) + '\n'))
    f.write(str(''.join(comment2) + '\n'))
    f.write(str(''.join(comment3) + '\n'))
    f.write(str(''.join(comment4) + '\n'))
    df_mode.to_csv(f, sep='\t', header=None, index=False)
    f.close()

    #Convert SAM to BAM and index
    os.system('mv ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.txt ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.sam')
    os.system('samtools view -S -b ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.sam > ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.bam')
    os.system('rm ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.sam')
    os.system('samtools sort ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode.bam -o ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode_sorted.bam')
    os.system('samtools index ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode_sorted.bam')

    #Perform plastid analysis
    os.system('metagene count -q ' + str(args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]) + 'metagene_reference_rois.txt ' + str(args_dict['metrics']) + str(file[:-4]) + '_periodicity --count_files ' + str(args_dict['metrics']) + str(file[:-4]) + '_mode_sorted.bam --fiveprime --offset 14 --normalize_over 30 200 --min_counts 50 --cmap Blues --title ' + file[:-4])

def make_periodicity(args_dict):

    args_dict = add_directory(args_dict, 'output', 'periodicity')
    args_dict = add_directory(args_dict, 'periodicity', 'metrics')
    output_location = args_dict['periodicity']

    #Get list of sam files from user input
    files = get_files(args_dict['input'], ['.sam'])

    #Add metagene reference to folder
    if 'generate_ref' in args_dict and args_dict['generate_ref'] == True:
        os.system('metagene generate ' + str(args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]) + 'metagene_reference --landmark ' + str(args_dict['landmark']) + ' --annotation_files ' + str(args_dict['gtf']) + ' --downstream ' + str(args_dict['downstream']))
    else:
        print('****************************************\nSkipping periodicity reference generation...\nIf parameters other than default values for landmark and downstream were used during reference generation, these values need to be provided even if not creating a reference...\n****************************************')

    #Perform metagene analysis
    parallelize(get_peaks, files, args_dict)

    #Compile images
    files = get_files(args_dict['metrics'], ['_periodicity_metagene_profile.txt'])
    compile_size_distribution(args_dict, args_dict['metrics'], files, 'metagene_average', str(int(args_dict['downstream'] - 1)), 'position', 'metagene coverage', 'periodicity', args_dict['experiment'], output_location)

    #Clean output
    os.system('rm ' + str(args_dict['metrics']) + '*_mode*')
    os.system('rm ' + str(args_dict['metrics']) + '*png')

"""
DESCRIPTION: Generate metagene profiles
"""
def get_profiles(args):

    file, args_dict = args[0], args[1]

    os.system('picard CollectRnaSeqMetrics REF_FLAT=' + str(args_dict['flat_type']) + ' STRAND_SPECIFICITY=NONE INPUT=' + str(args_dict['input']) + str(file) + ' OUTPUT=' + str(args_dict['metrics']) + str(file[:-4]) + '_rna_metrics')

def make_metagene(args_dict):

    #Add output directory to output for metagene profiles
    args_dict = add_directory(args_dict, 'output', 'metagene')
    args_dict = add_directory(args_dict, 'metagene', 'metrics')
    output_location = args_dict['metagene']

    #Get list of sam files from user input
    files = get_files(args_dict['input'], ['.sam'])

    #Perform metagene analysis
    parallelize(get_profiles, files, args_dict)

    #Compile images
    files = get_files(args_dict['metrics'], ['_rna_metrics'])
    compile_size_distribution(args_dict, args_dict['metrics'], files, 'normalized_position', None, 'meta_position', 'normalized_coverage (all_reads)', 'metagene',  args_dict['experiment'], output_location)

"""
DESCRIPTION: Generate read distribution profiles
"""
def run_fastqc(args):

    file, args_dict = args[0], args[1]

    os.system("fastqc -q " + str(args_dict['input']) + str(file) + " -o " + str(args_dict['fastqc_output']))

def make_readDistributions(args_dict):

    args_dict = add_directory(args_dict, 'output', 'read_distributions')
    args_dict = add_directory(args_dict, 'read_distributions', 'fastqc_output')
    output_location = args_dict['read_distributions']

    #Get FASTQC file list and unzip
    files = get_files(args_dict['input'], ['.fastq', '.fq', '.txt'])

    #Perform fastqc on each file and unzip output
    parallelize(run_fastqc, files, args_dict)

    files = get_files(args_dict['fastqc_output'], ['.zip'])
    for file in files:
        if file.endswith('.zip'):
            os.system('unzip -n -q ' + str(args_dict['fastqc_output']) + str(file) + ' -d ' + str(args_dict['fastqc_output']))

    #Compile read distributions
    files = get_files(args_dict['fastqc_output'], ['.zip'])
    files = [str(x[:-4]) + '/fastqc_data.txt' for x in files]
    compile_size_distribution(args_dict, args_dict['fastqc_output'], files, '#Length', '>>END_MODULE', 'position', 'read size (bp)', 'fastqc', args_dict['experiment'], output_location)
