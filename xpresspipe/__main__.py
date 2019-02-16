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
from .__init__ import __version__
from .messages import *
from .arguments import get_arguments
from .trim import run_trim
from .align import run_seRNAseq, run_peRNAseq
from .count import create_bed, create_bigwig, count_reads, collect_counts
from .quality import *
from .reference import truncate, rrnaprobe

"""
DESCRIPTION: Main function to call necessary functions for sub-modules

ASSUMPTIONS:
Proper arguments are provided where some user renaming of files may be required
"""
def main(args=None):

    #Print license information
    msg_license()

    #Collate CLI arguments provided by user
    try:
        args, args_dict = get_arguments(args, __version__)
    except:
        raise Exception("There was an issue in processing the arguments for the provided function.")

    #Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'trim':
        run_trim(args_dict)

    elif args.cmd == 'align':
        #Align
        if args_dict['type'].upper() == 'SE':
            args_dict = run_seRNAseq(args_dict)
        elif args_dict['type'].upper() == 'PE':
            args_dict = run_peRNAseq(args_dict)
        else:
            raise Exception('Invalid type argument provided')
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])

    elif args.cmd == 'count':
        print('coming soon')

    elif args.cmd == 'quality':
        print('coming soon')

    elif args.cmd == 'truncate':
        print('coming soon')

    elif args.cmd == 'rrnaProbe':
        print('coming soon')

    elif args.cmd == 'convertNames':
        print('coming soon')

    elif args.cmd == 'normalizeTable':
        print('coming soon')

    elif args.cmd == 'seRNAseq':
        #Trim
        args_dict = run_trim(args_dict)
        #Align
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count
        args_dict['input'] = args_dict['aligndir']


    elif args.cmd == 'peRNAseq':
        #Trim
        args_dict = run_trim(args_dict)
        #Align
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_peRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count
        args_dict['input'] = args_dict['aligndir']


    elif args.cmd == 'riboprof':
        #Trim
        args_dict = run_trim(args_dict)
        #Align
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count
        args_dict['input'] = args_dict['aligndir']


    else:
        raise Exception("Invalid function processing function provided.")

"""
DESCRIPTION: Run main
"""
if __name__ == "__main__":

    sys.exit(main() or 0)
