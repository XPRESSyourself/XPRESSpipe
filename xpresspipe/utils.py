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

"""Check directory formatting"""
def check_directories(input):

    # Check input directory name is formatted correctly and fix if necessary
    if input.endswith('/'):
        pass
    elif '.' in input: # Assumes this is a file input and not a directory input
        pass
    else:
        input += '/'

    return input

"""Create output directory"""
def add_directory(
    args_dict, parent, name):

    os.system('mkdir' \
        + ' ' + str(args_dict[str(parent)]) + str(name) \
        + str(args_dict['log']))

    args_dict[name] = str(str(args_dict[str(parent)]) + str(name) + '/')

    return args_dict

"""Make a list of the files in a given directory, based on list of acceptable file suffixes"""
def get_files(
    directory, suffix, omit=[]):

    # Initialize blank file list to fill
    file_list = []

    # Walk through raw data files within given directory
    for file in os.listdir(directory):

        for s in suffix:
            if file.endswith(str(s)):
                file_list.append(file) # Do not append directory, files in list will be modified and output to different locations

    # Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    # Get rid of bad grabs
    omit_drop = []
    if len(omit) > 0:

        for x in file_list:

            for o in omit:
                if str(o) in x:
                    omit_drop.append(x)

    for x in omit_drop:
        file_list.remove(x)

    return tuple(file_list)

"""Get files to perform rRNA prober upon"""
def get_probe_files(
    args_dict, suffix):

    # Initialize blank file list to fill
    probe_list = []

    # Walk through raw data files within given directory
    for file in os.listdir(args_dict['input']):

        for s in suffix:
            if file.endswith(str(s)):
                probe_list.append(file)

    for x in args_dict['input']:
        if x.endswith(str(suffix)) == True:
            if 'footprint_only' in args_dict:
                if 'FOOTPRINT' in x.upper() or 'FP' in x.upper() or 'RPF' in x.upper():
                    probe_list.append(args_dict['input'] + x)
            else:
                probe_list.append(args_dict['input'] + x)

    print(probe_list)

    return tuple(probe_list)

"""Unzip all files from directory"""
def unzip_files(directory):

    suffix = ['.gz', '.zip']

    # Walk through raw data files within given directory
    for file in os.listdir(directory):

        for s in suffix:
            if file.endswith(str(s)):
                if s == '.gz':
                    os.system('gzip' +
                        ' -d' \
                        + ' ' + str(directory) + str(file) \
                        + str(args_dict['log']))
                if s == '.zip':
                    os.system('unzip' \
                        + ' ' + str(directory) + str(file) \
                        + str(args_dict['log']))

"""Get fasta files within directory"""
def get_fasta(fasta_directory):

    # Make space separated list of fasta files
    fasta = get_files(fasta_directory, ['.txt', '.fasta', '.fa'], omit=['refFlat', 'rois'])
    fasta = [fasta_directory + x for x in fasta]

    if len(fasta) > 1:
        fasta_list = ' '.join(fasta)
    else:
        fasta_list = ''.join(fasta)

    return fasta_list

"""Make periodicity flat reference using riboWaltz"""
def make_flat(args_dict):

    if not str(args_dict['gtf']).endswith('.gtf'):
        raise Exception('GTF file required')

    os.system('rscript' \
        + ' ' + str(__path__) + '/periodicity.r' \
        + ' TRUE' \
        + ' ' + str(args_dict['gtf']) \
        + ' NONE' \
        + ' NONE' \
        + str(args_dict['log']))
