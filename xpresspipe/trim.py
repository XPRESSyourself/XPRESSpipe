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
from .utils import get_files
from .parallel import parallelize, parallelize_pe

"""
DESCRIPTION: Add output directories
"""
def add_trim_directories(args_dict):

    os.system('mkdir ' + str(args_dict['output']) + 'trimmed_fastq')
    args_dict['trimdir']: str(str(args_dict['output']) + 'trimmed_fastq/'),

    return args_dict

"""
DESCRIPTION: Determine sequencing type based on adaptor list
"""
def determine_type(adaptor_list):

    try:
        if adaptor_list == None:
            return None
        if len(adaptor_list) == 1:
            if adaptor_list[0].upper() == 'POLYX':
                return 'POLYX'
            else:
                return 'SE'
        if len(adaptor_list) == 2:
            return 'PE'

    except:
        raise Exception('Could not determine sequencing type from adaptor list during read trimming')

"""
DESCRIPTION: Trim adaptors via auto-detect
"""
def auto_trim(args):

    file, args_dict = args[0], args[1]
    os.system("fastp -f 1 -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "trimmed_" + file + " -l " + str(args_dict['read_length_min']) + " -q " + str(args_dict['read_quality']) + " -j " + str(args_dict['trimdir']) + file[:-6] + "fastp.json -h " + str(args_dict['trimdir']) + file[:-6] + "fastp.html")

"""
DESCRIPTION: Trim polyX adaptors
"""
def polyx_trim(args):

    file, args_dict = args[0], args[1]
    os.system("fastp -f 1 -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "trimmed_" + file + " --polyX -l " + str(args_dict['read_length_min']) + " -q " + str(args_dict['read_quality']) + " -j " + str(args_dict['trimdir']) + file[:-6] + "fastp.json -h " + str(args_dict['trimdir']) + file[:-6] + "fastp.html")

"""
DESCRIPTION: Trim SE adaptors
"""
def se_trim(args):

    file, args_dict = args[0], args[1]
    os.system("fastp -f 1 -i " + str(args_dict['input']) + file + " -o " + str(args_dict['output']) + "trimmed_" + file + " -a " + str(args_dict['adaptors'][0]) + " -l " + str(args_dict['read_length_min']) + " -q " + str(args_dict['read_quality']) + " -j " + str(args_dict['trimdir']) + file[:-6] + "fastp.json -h " + str(args_dict['trimdir']) + file[:-6] + "fastp.html")

"""
DESCRIPTION: Trim PE adaptors
"""
def pe_trim(args):

    file1, file2, args_dict = args[0], args[1], args[2]
    os.system("fastp -f 1 -i " + str(args_dict['input']) + file1 + " -I " + str(args_dict['input']) + file2 + " -o " + str(args_dict['output']) + "trimmed_" + file1 + " -O " + str(args_dict['output']) + "trimmed_" + file2 + " -a " + str(args_dict['adaptors'][0]) + " --adapter_sequence_r2 " + str(args_dict['adaptors'][1]) + " -l " + str(args_dict['read_length_min']) + " -q " + str(args_dict['read_quality']) + " -j " + str(args_dict['trimdir']) + file[:-6] + "fastp.json -h " + str(args_dict['trimdir']) + file[:-6] + "fastp.html")

"""
DESCRIPTION: Trim RNAseq reads of adaptors and for quality
"""
def run_trim(args_dict):

    try:
        #Add output directories
        args_dict = add_trim_directories(args_dict)

        #Get list of files to trim based on acceptable file types
        files = get_files(args_dict['input'], ['.fastq','.fq','.txt'])

        #Determine sequencing type based on args_dict['adaptors']
        type = determine_type(args_dict['adaptors'])

        #Auto-detect and trim adaptors
        if type == None:
            parallize(auto_trim, files, args_dict)

        #Trim polyX adaptors, assumes single-end RNAseq reads
        if type == 'POLYX':
            parallize(polyx_trim, files, args_dict)

        #Trim single-end RNAseq reads
        if type == 'SE':
            parallize(se_trim, files, args_dict)

        #Trim paired-end RNAseq reads
            parallelize_pe(pe_trim, files, args_dict)

    except:
        raise Exception('Trimming failed')
