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

"""
DESCRIPTION: Check directory formatting
"""
#Check directory formatting
def check_directories(directory):

    #Check input directory name is formatted correctly and fix if necessary
    if directory.endswith('/'):
        pass
    else:
        directory += '/'

    return directory

"""
DESCRIPTION: Make a list of the files in a given directory, based on list of acceptable file suffixes
"""
def get_files(directory, suffix):

    #Initialize blank file list to fill
    file_list = []

    #Walk through raw data files within given directory
    for subdir, dirs, files in os.walk(directory):
        for f in files:
            for s in suffix:
                if f.endswith(str(s)):
                    file_list.append(f)

    #Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    return tuple(file_list)

"""
"""
def get_probe_files(args_dict, suffix):

    #Initialize blank file list to fill
    probe_list = []

    #Walk through raw data files within given directory
    for subdir, dirs, files in os.walk(args_dict['input']):
        for f in files:
            for s in suffix:
                if f.endswith(str(s)):
                    file_list.append(f)


    for x in args_dict['input']:
        if x.endswith(str(suffix)) == True:
            if 'footprint_only' in args_dict:
                if 'FOOTPRINT' in x.upper() or 'FP' in x.upper() or 'RPF' in x.upper():
                    probe_list.append(args_dict['input'] + x)
            else:
                probe_list.append(args_dict['input'] + x)

    return tuple(probe_list)

"""
DESCRIPTION: Unzip all files from directory
"""
def unzip_files(directory):

    suffix = ['.gz', '.zip']

    #Walk through raw data files within given directory
    for subdir, dirs, files in os.walk(directory):
        for f in files:
            for s in suffix:
                if f.endswith(str(s)):
                    if s == '.gz':
                        os.system('gzip -d ' + str(directory) + str(f))
                    if s == '.zip':
                        os.system('unzip ' + str(directory) + str(f))
