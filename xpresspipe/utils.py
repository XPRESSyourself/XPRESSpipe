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

"""Check directory formatting"""
def check_directories(
    input,
    type=None):

    # Check that a file wasn't passed in
    if os.path.isdir(input) != True:
        raise Exception(str(input) + ' is not a directory')

    # Check input directory name is formatted correctly and fix if necessary
    input = os.path.abspath(input)

    if input.endswith('/'):
        pass
    else:
        input += '/'

    return input

"""Create output directory"""
def add_directory(
    args_dict,
    parent,
    name):

    # Check that name is valid
    if '.' in name:
        raise Exception('Invalid new directory name provided')

    if name.endswith('/'):
        name = name[:-1]

    os.system(
        'mkdir'
        + ' -p'
        + ' ' + str(args_dict[str(parent)]) + str(name)
        + str(args_dict['log']))

    args_dict[name] = str(str(args_dict[str(parent)]) + str(name) + '/')

    return args_dict

"""Make a list of the files in a given directory, based on list of acceptable file suffixes"""
def get_files(
    directory,
    suffix,
    omit=[]):

    # Initialize blank file list to fill
    file_list = []

    # Walk through raw data files within given directory
    for file in next(os.walk(directory))[2]:

        if file.endswith(tuple(suffix)) and os.path.isfile(str(directory) + str(file)) == True:
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

    if len(file_list) == 0:
        print('The provided directory does not contain any files with the suffix ' + str(''.join(suffix)) + '.\nPlease check the suffixes of your files and see they are different than ' + str(''.join(suffix)) + '. If working with BAM files, try using the --bam_suffix argument to specify a different BAM suffix.')
        sys.exit(1)

    return tuple(file_list)

"""Make a list of the directories in a given directory, based on list of acceptable suffixes"""
def get_directories(
    directory,
    suffix,
    omit=[]):

    # Initialize blank file list to fill
    directory_list = []

    # Walk through raw data files within given directory
    for dir in next(os.walk(directory))[1]:
        if dir.endswith(tuple(suffix)):
            directory_list.append(dir) # Do not append directory, files in list will be modified and output to different locations

    # Get full path for each directory
    directory_list = [str(directory) + str(x) + '/' for x in directory_list]

    # Sort files in alphabetical order (helps in formatting the count tables correctly)
    directory_list = sorted(directory_list)

    # Get rid of bad grabs
    omit_drop = []
    if len(omit) > 0:
        for x in directory_list:
            for o in omit:
                if str(o) in x:
                    omit_drop.append(x)

    for x in omit_drop:
        directory_list.remove(x)

    return tuple(directory_list)

"""Get files to perform rRNA prober upon"""
"""def get_probe_files(
    args_dict,
    suffix):

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

    return tuple(probe_list)"""

"""Unzip all files from directory"""
def unzip_files(
    directory,
    log):

    suffix = ['.gz', '.zip']

    # Walk through raw data files within given directory
    for file in os.listdir(directory):

        for s in suffix:
            if file.endswith(str(s)):
                if s == '.gz':
                    os.system(
                        'gzip' +
                        ' -d'
                        + ' ' + str(directory) + str(file)
                        + str(log))

                if s == '.zip':
                    os.system(
                        'unzip'
                        + ' -q'
                        + ' -d ' + str(directory)
                        + ' ' + str(directory) + str(file)
                        + str(log))

"""Check path for directory named fasta with fasta files"""
"""
- A one-level recursive function for searching parent and child directories for fasta files
- Get list of parent directory
- If a list of files that exist cannot be found, recurse through child directories until list found
    - If a list can be found, return list
    - If list cannot be found in level-1 children directories, exit and print Error
"""
def get_fasta(
    fasta_directory,
    suffix=['.fasta', '.fa'],
    omit=['refFlat', 'rois'],
    recursed=0):

    fasta_list = []
    fasta_directory = check_directories(fasta_directory)

    # Walk through raw data files within given directory
    for file in next(os.walk(fasta_directory))[2]:
        if file.endswith(tuple(suffix)) and os.path.isfile(str(fasta_directory) + str(file)) == True:
            fasta_list.append(file) # Do not append directory, files in list will be modified and output to different locations

    # Sort files in alphabetical order (helps in formatting the count tables correctly)
    fasta_list = sorted(fasta_list)

    # Get rid of bad grabs
    omit_drop = []
    if len(omit) > 0:
        for x in fasta_list:
            for o in omit:
                if str(o) in x:
                    omit_drop.append(x)

    for x in omit_drop:
        fasta_list.remove(x)

    # If fasta list is empty, try to find a child directory with fasta files
    if len(fasta_list) == 0:
        print('FASTA files not found in parent reference directory ' + str(fasta_directory) + '.\nAttempting search through any child directories...')

        # Go through first level child directories until a valid list of fasta files is found
        while recursed == 0:
            for folder in os.listdir(fasta_directory):
                folder = str(fasta_directory) + str(folder)
                if os.path.isdir(folder) == True:
                    fasta_list = get_fasta(
                        folder,
                        recursed=recursed)

                    # Break out of the loop if valid fasta list found
                    if len(fasta_list) > 0:
                        break
                    else:
                        continue

            # Prep output -- return list if valid, or send output trigger for a sys.exit()
            if len(fasta_list) > 0:
                print('A list of valid FASTA files found in ' + str(folder))
                recursed = 2
            else:
                recursed = 1 # Prevent excessive recursions

        # What to do at the end -- either return the list or print exit message
        else:
            if recursed == 1:
                # No child directories have fasta files if this is hit
                print('A directory with fasta files cannot be found.\nExiting...')
                sys.exit(1)
            else:
                return fasta_list

    else:
        # If valid list of fastas available, create a single string out of the fullpath/filename of each fasta
        fasta = [fasta_directory + x for x in fasta_list]
        if len(fasta) > 1:
            fasta_list = ' '.join(fasta)
        else:
            fasta_list = ''.join(fasta)

        return fasta_list
