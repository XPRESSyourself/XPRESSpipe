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
from .utils import get_probe_files, create_reference, get_summary, create_flat
from .quality import make_metagene, make_readDistributions, make_periodicity
from xpresstools import truncate, rpm, r_fpkm, log_scale, batch_normalize, convert_names_gtf

"""
DESCRIPTION: Main function to call necessary functions for sub-modules

ASSUMPTIONS:
Proper arguments are provided where some user renaming of files may be required
"""
def main(args=None):

    #Print license information
    #msg_license()
    args, args_dict = get_arguments(args, __version__)

    #Get GTF type
    if 'count_coding' in args_dict:
        if 'truncate' in args_dict:
            args_dict['gtf_type'] = str(args_dict['reference']) + 'transcripts_coding_truncated.gtf'
            args_dict['flat_type'] = str(args_dict['reference']) + 'transcripts_coding_truncated_refFlat.txt'
        else:
            args_dict['gtf_type'] = str(args_dict['reference']) + 'transcripts_coding.gtf'
            args_dict['flat_type'] = str(args_dict['reference']) + 'transcripts_coding_refFlat.txt'
    else:
        args_dict['gtf_type'] = str(args_dict['reference']) + 'transcripts.gtf'
        args_dict['flat_type'] = str(args_dict['reference']) + 'transcripts_refFlat.txt'

    #Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'trim':
        print('Trimming reads...')
        run_trim(args_dict)
        msg_complete()

    elif args.cmd == 'align':
        #Align
        print('Aligning reads to reference...')
        if args_dict['type'].upper() == 'SE':
            args_dict = run_seRNAseq(args_dict)
        elif args_dict['type'].upper() == 'PE':
            args_dict = run_peRNAseq(args_dict)
        else:
            raise Exception('Invalid type argument provided')
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['alignments'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['alignments'])
        msg_complete()

    elif args.cmd == 'count':
        print('Counting alignments...')
        #Count reads for each alignment file
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        print('Collecting and collating counts...')
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        msg_complete()

    elif args.cmd == 'metagene':
        print('Performing metagene analysis on SAM files...')
        make_metagene(args_dict)
        msg_complete()

    elif args.cmd == 'readDistribution':
        print('Performing read distribution analysis on fastq files...')
        make_readDistributions(args_dict)
        msg_complete()

    elif args.cmd == 'periodicity':
        print('Performing periodicity analysis on most abundant read length in SAM files...')
        make_periodicity(args_dict)
        msg_complete()

    elif args.cmd == 'truncate':
        print('Formatting coding only and truncated reference files...')
        output_path = args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]
        truncate(args_dict['gtf'], truncate_amount=args_dict['truncate_amount'], save_coding_path=str(output_path), save_truncated_path=str(output_path), sep='\t', return_files=False)
        msg_complete()

    elif args.cmd == 'makeFlat':
        print('Formatting coding only and truncated reference files...')
        create_flat(args_dict['input'])
        msg_complete()

    elif args.cmd == 'createReference':
        print('Creating reference files...')
        create_reference(args_dict['output'], args_dict['fasta'], args_dict['gtf'], threads=args_dict['threads'], sjdbOverhang=args_dict['sjdbOverhang'])
        msg_complete()

    elif args.cmd == 'rrnaProbe':
        #Get files to probe
        probe_list = get_probe_files(args_dict, '.zip')
        #Run rrna_prober, output to outputDir
        print('Probing for most over-represented read sequences...')
        probe_out = rrnaProbe(probe_list, args_dict['min_overlap']) #use inputDir to get FASTQC files and output to outputDir/analysis
        #Output summary
        with open(args_dict['output'] + 'rrnaProbe_output.txt', "w") as text_file:
            print(probe_out, file=text_file)
        msg_complete()

    elif args.cmd == 'convertNames':
        #Convert row names in dataframe
        print('Converting row names...')
        convert_names_gtf(args_dict['data'], args_dict['gtf'], orig_name_label=args_dict['orig_name_label'], orig_name_location=args_dict['orig_name_location'], new_name_label=args_dict['new_name_label'], new_name_location=args_dict['new_name_location'], refill=args_dict['refill'], sep='\t')
        msg_complete()

    elif args.cmd == 'normalizeMatrix':
        #Run in sample normalization
        print('Normalizing matrix...')
        run_normalization(args_dict)
        msg_complete()

    elif args.cmd == 'seRNAseq':
        args_dict['type'] = 'SE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['alignments'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['alignments'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['alignments']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        #Normalize
        msg_normalize()
        args_dict['data'] = str(args_dict['experiment']) + 'counts_table.csv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_metagene(args_dict)

        get_summary(args_dict)
        msg_finish()

    elif args.cmd == 'peRNAseq':
        args_dict['type'] = 'PE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_peRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['alignments'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['alignments'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['alignments']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        #Normalize
        msg_normalize()
        args_dict['data'] = str(args_dict['experiment']) + 'counts_table.csv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_metagene(args_dict)

        get_summary(args_dict)
        msg_finish()

    elif args.cmd == 'riboprof':
        args_dict['type'] = 'SE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        if args_dict['output_bed'] == True:
            create_bed(args_dict['output'], args_dict['alignments'])
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict['output'], args_dict['alignments'])
        #Count reads for each alignment file
        msg_count()
        args_dict['input'] = args_dict['alignments']
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        #Normalize
        msg_normalize()
        args_dict['data'] = str(args_dict['experiment']) + 'counts_table.csv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_periodicity(args_dict)
        make_metagene(args_dict)

        get_summary(args_dict)
        msg_finish()

    else:
        raise Exception("Invalid function processing function provided.")

"""
DESCRIPTION: Run main
"""
if __name__ == "__main__":

    sys.exit(main() or 0)
