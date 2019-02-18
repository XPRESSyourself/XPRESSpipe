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
from .count import create_bed, create_bigwig, count_reads, collect_counts, run_normalization
from .rrnaprobe import rrnaProbe
from .utils import get_probe_files
from .quality import *#FILL THESE IN ONCE CREATED
from xpresstools import truncate, rpm, r_fpkm, log_scale, batch_normalize, convert_names_gtf

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

    #Sort files in

    #Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'trim':
        msg_trim()
        run_trim(args_dict)

    elif args.cmd == 'align':
        #Align
        msg_align()
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
        msg_count()
        #Count reads for each alignment file
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)

    elif args.cmd == 'quality':
        print('coming soon')

    elif args.cmd == 'truncate':
        output_path = args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]
        truncate(args_dict['gtf'], truncate_amount=args_dict['truncate_amount'], save_coding_path=str(output_path), save_truncated_path=str(output_path), sep='\t', return_files=False)

    elif args.cmd == 'rrnaProbe':
        #Get files to probe
        probe_list = get_probe_files(args_dict, '.zip')
        #Run rrna_prober, output to outputDir
        probe_out = rrnaProbe(probe_list, args_dict['min_overlap']) #use inputDir to get FASTQC files and output to outputDir/analysis
        #Output summary
        with open(args_dict['output'] + 'rrnaProbe_output.txt', "w") as text_file:
            print(probe_out, file=text_file)

    elif args.cmd == 'convertNames':
        #Convert row names in dataframe
        convert_names_gtf(args_dict['data'], args_dict['gtf'], orig_name_label=args_dict['orig_name_label'], orig_name_location=args_dict['orig_name_location'], new_name_label=args_dict['new_name_label'], new_name_location=args_dict['new_name_location'], refill=args_dict['refill'], sep='\t')

    elif args.cmd == 'normalizeMatrix':
        #Run in sample normalization
        run_normalization(args_dict)

    elif args.cmd == 'seRNAseq':
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['aligndir']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        matrix_location = str(args_dict['experiment']) + 'counts_table.csv'
        #Normalize
        args_dict['data'] = matrix_location
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        run_normalization(args_dict)

        msg_finish()

    elif args.cmd == 'peRNAseq':
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_peRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['aligndir']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        matrix_location = str(args_dict['experiment']) + 'counts_table.csv'
        #Normalize
        args_dict['data'] = matrix_location
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        run_normalization(args_dict)

        msg_finish()

    elif args.cmd == 'riboprof':
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimdir']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['aligndir'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['aligndir'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['aligndir']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        matrix_location = str(args_dict['experiment']) + 'counts_table.csv'
        #Normalize
        args_dict['data'] = matrix_location
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        run_normalization(args_dict)

        msg_finish()

    else:
        raise Exception("Invalid function processing function provided.")

"""
DESCRIPTION: Run main
"""
if __name__ == "__main__":

    sys.exit(main() or 0)
