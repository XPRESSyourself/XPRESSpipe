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

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, add_directory
from .parallel import parallelize, parallelize_pe

"""Determine sequencing type based on adaptor list"""
def determine_type(
        adaptor_list):

    try:
        if adaptor_list == None:
            return 'AUTOSE'
        else:
            # Convert case
            adaptor_list = [x.upper() for x in adaptor_list]
            if len(adaptor_list) == 1:
                if 'NONE' in adaptor_list:
                    return 'AUTOSE'
                elif 'POLYX' in adaptor_list:
                    return 'POLYX'
                else:
                    return 'SE'
            elif len(adaptor_list) == 2:
                if 'NONE' in adaptor_list:
                    return 'AUTOPE'
                else:
                    return 'PE'
            else:
                raise Exception('Invalid number of adaptor options')

    except:
        raise Exception('Could not determine sequencing type from adaptor list during read trimming')

"""Trim adaptors via auto-detect"""
def auto_trim(
        args):

    file, args_dict = args[0], args[1] # Parse args

    os.system(
        'fastp'
        + ' -f 1'
        + ' --thread ' + str(args_dict['threads'])
        + ' -i ' + str(args_dict['input']) + str(file)
        + ' -o ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file)
        + ' -l ' + str(args_dict['min_length'])
        + ' -q ' + str(args_dict['quality'])
        + ' -j ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.json'
        + ' -h ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.html'
        + str(args_dict['log']))

"""Trim polyX adaptors"""
def polyx_trim(
        args):

    file, args_dict = args[0], args[1] # Parse args

    os.system(
        'fastp'
        + ' -f 1'
        + ' --thread ' + str(args_dict['threads'])
        + ' -i ' + str(args_dict['input']) + file
        + ' -o ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file)
        + ' --trim_poly_x'
        + ' -l ' + str(args_dict['min_length'])
        + ' -q ' + str(args_dict['quality'])
        + ' -j ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.json'
        + ' -h ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.html'
        + str(args_dict['log']))

"""Trim SE adaptors"""
def se_trim(
        args):

    file, args_dict = args[0], args[1] # Parse args

    os.system(
    'fastp'
    + ' -f 1'
    + ' --thread ' + str(args_dict['threads'])
    + ' -i ' + str(args_dict['input']) + str(file)
    + ' -o ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file)
    + ' -a ' + str(args_dict['adaptors'][0])
    + ' -l ' + str(args_dict['min_length'])
    + ' -q ' + str(args_dict['quality'])
    + ' -j ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.json'
    + ' -h ' + str(args_dict['trimmed_fastq']) + str(file[:-6]) + 'fastp.html'
    + str(args_dict['log']))

"""Auto-detect and trim PE adaptors"""
def auto_pe_trim(
        args):

    file1, file2, args_dict = args[0], args[1], args[2] # Parse args

    os.system(
        'fastp'
        + ' -f 1'
        + ' --thread ' + str(args_dict['threads'])
        + ' -i ' + str(args_dict['input']) + str(file1)
        + ' -I ' + str(args_dict['input']) + str(file2)
        + ' -o ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file1)
        + ' -O ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file2)
        + ' -l ' + str(args_dict['min_length'])
        + ' -q ' + str(args_dict['quality'])
        + ' -j ' + str(args_dict['trimmed_fastq']) + str(file1[:-6]) + 'fastp.json'
        + ' -h ' + str(args_dict['trimmed_fastq']) + str(file1[:-6]) + 'fastp.html'
        + str(args_dict['log']))

"""Trim PE adaptors"""
def pe_trim(
        args):

    file1, file2, args_dict = args[0], args[1], args[2] # Parse args

    os.system(
        'fastp'
        + ' -f 1'
        + ' --thread ' + str(args_dict['threads'])
        + ' -i ' + str(args_dict['input']) + str(file1)
        + ' -I ' + str(args_dict['input']) + str(file2)
        + ' -o ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file1)
        + ' -O ' + str(args_dict['trimmed_fastq']) + 'trimmed_' + str(file2)
        + ' -a ' + str(args_dict['adaptors'][0]) + ' --adapter_sequence_r2 ' + str(args_dict['adaptors'][1])
        + ' -l ' + str(args_dict['min_length'])
        + ' -q ' + str(args_dict['quality'])
        + ' -j ' + str(args_dict['trimmed_fastq']) + str(file1[:-6]) + 'fastp.json'
        + ' -h ' + str(args_dict['trimmed_fastq']) + str(file1[:-6]) + 'fastp.html'
        + str(args_dict['log']))

"""Trim RNAseq reads of adaptors and for quality"""
def run_trim(
        args_dict):

    try:
        # Add output directories
        args_dict = add_directory(
            args_dict,
            'output',
            'trimmed_fastq')

        # Get list of files to trim based on acceptable file types
        files = get_files(
            args_dict['input'],
            ['.fastq','.fq','.txt'])

        # Determine sequencing type based on args_dict['adaptors']
        if 'type' not in args_dict:
            type = determine_type(
                args_dict['adaptors'])
        else:
            type = args_dict['type']

        # Auto-detect and trim adaptors
        if type == 'AUTOSE':
            parallelize(
                auto_trim,
                files,
                args_dict)

        # Trim polyX adaptors, assumes single-end RNAseq reads
        if type == 'POLYX':
            parallelize(
                polyx_trim,
                files,
                args_dict)

        # Trim single-end RNAseq reads
        if type == 'SE':
            parallelize(
                se_trim,
                files,
                args_dict)

        # Auto-detect adaptors for paired-end reads
        if type == 'AUTOPE':
            if len(files) % 2 != 0:
                raise Exception('An uneven number of paired-end files were specified in the input directory')
            else:
                parallelize_pe(
                    auto_pe_trim,
                    files,
                    args_dict)

        # Trim paired-end RNAseq reads
        if type == 'PE':
            if len(files) % 2 != 0:
                raise Exception('An uneven number of paired-end files were specified in the input directory')
            else:
                parallelize_pe(
                    pe_trim,
                    files,
                    args_dict)

        return args_dict

    except:
        raise Exception('Trimming failed')
