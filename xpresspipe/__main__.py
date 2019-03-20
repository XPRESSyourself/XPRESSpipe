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
import pandas as pd
from .__init__ import __version__
from .messages import *
from .arguments import get_arguments
from .trim import run_trim
from .align import run_seRNAseq, run_peRNAseq
from .count import create_bed, create_bigwig, count_reads, collect_counts, run_normalization
from .rrnaprobe import rrnaProbe
from .utils import get_probe_files, create_reference, get_summary, create_flat, unzip_files
from .quality import make_metagene, make_readDistributions, make_periodicity
from .parallel import get_cores
from xpresstools import truncate, rpm, r_fpkm, log_scale, batch_normalize, convert_names_gtf, diff_xpress

"""
DESCRIPTION: Main function to call necessary functions for sub-modules

ASSUMPTIONS:
Proper arguments are provided where some user renaming of files may be required
"""
def main(args=None):

    #Print license information
    #msg_license()
    args, args_dict = get_arguments(args, __version__)

    #Should have already seen check_directory() so should have a trailing '/'
    if 'input' in args_dict and str(args_dict['input']).endswith('/'):
        unzip_files(args_dict['input'])

    #Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'trim':
        print('Trimming reads...')
        run_trim(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'TRIM')

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
        args_dict['input'] = args_dict['alignments']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'ALIGN')

    elif args.cmd == 'count':
        print('Counting alignments...')
        #Count reads for each alignment file
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        print('Collecting and collating counts...')
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'COUNT')

    elif args.cmd == 'diffxpress':
        print('Performing differential expression analysis...')
        diff_xpress(str(args_dict['input']), str(args_dict['sample']), equation=str(args_dict['design']))
        check_process(args_dict['log_file'], msg_complete(), 'DIFFERENTIAL EXPRESSION')

    elif args.cmd == 'metagene':
        print('Performing metagene analysis on SAM files...')
        make_metagene(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'METAGENE')

    elif args.cmd == 'readDistribution':
        print('Performing read distribution analysis on fastq files...')
        make_readDistributions(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'READ DISTRIBUTION')

    elif args.cmd == 'periodicity':
        print('Performing periodicity analysis on most abundant read length in SAM files...')
        make_periodicity(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'PERIODICITY')

    elif args.cmd == 'curateReference':
        print('Curating reference')
        args_dict['create_refFlats'] = True
        #Create STAR reference
        args_dict['threads'], args_dict['workers'] = get_cores(args_dict)
        create_reference(args_dict['output'], args_dict['fasta'], args_dict['gtf'], threads=args_dict['threads'], sjdbOverhang=args_dict['sjdbOverhang'])
        #Truncate transcript reference
        truncate(args_dict['gtf'], truncate_amount=args_dict['truncate_amount'], save_coding_path=str(args_dict['output']), save_truncated_path=str(args_dict['output']), sep='\t', return_files=False)
        #Flatten transcript reference files
        create_flat(args_dict['output'])
        check_process(args_dict['log_file'], msg_complete(), 'CURATE REFERENCE')

    elif args.cmd == 'makeReference':
        print('Creating reference files...')
        args_dict['threads'], args_dict['workers'] = get_cores(args_dict)
        create_reference(args_dict['output'], args_dict['fasta'], args_dict['gtf'], threads=args_dict['threads'], sjdbOverhang=args_dict['sjdbOverhang'])
        check_process(args_dict['log_file'], msg_complete(), 'MAKE REFERENCE')

    elif args.cmd == 'truncate':
        print('Formatting coding only and truncated reference files...')
        output_path = args_dict['gtf'][:args_dict['gtf'].rfind('/') + 1]
        truncate(args_dict['gtf'], truncate_amount=args_dict['truncate_amount'], save_coding_path=str(output_path), save_truncated_path=str(output_path), sep='\t', return_files=False)
        if 'create_refFlats' in args_dict and args_dict['create_refFlats'] == True:
            create_flat(str(output_path))
        check_process(args_dict['log_file'], msg_complete(), 'TRUNCATE')

    elif args.cmd == 'makeFlat':
        print('Formatting coding only and truncated reference files...')
        create_flat(args_dict['input'])
        check_process(args_dict['log_file'], msg_complete(), 'MAKE REFFLAT')

    elif args.cmd == 'rrnaProbe':
        #Get files to probe
        probe_list = get_probe_files(args_dict, '.zip')
        #Run rrna_prober, output to outputDir
        print('Probing for most over-represented read sequences...')
        probe_out = rrnaProbe(probe_list, args_dict['min_overlap']) #use inputDir to get FASTQC files and output to outputDir/analysis
        #Output summary
        with open(args_dict['output'] + 'rrnaProbe_output.txt', "w") as text_file:
            print(probe_out, file=text_file)
        check_process(args_dict['log_file'], msg_complete(), 'RRNA PROBE')

    elif args.cmd == 'convertNames':
        #Convert row names in dataframe
        print('Converting row names...')
        if str(args_dict['input']).endswith('.csv'):
            delim = ','
            suf = '.csv'
        else:
            delim = '\t'
            suf = '.tsv'
        data = pd.read_csv(str(args_dict['input']), sep=delim)
        data = convert_names_gtf(data, args_dict['gtf'], orig_name_label=args_dict['orig_name_label'], orig_name_location=args_dict['orig_name_location'], new_name_label=args_dict['new_name_label'], new_name_location=args_dict['new_name_location'], refill=args_dict['refill'], sep='\t')
        data.to_csv(str(args_dict['input'])[:-4] + '_renamed' + str(suf), sep=delim, index=False)
        check_process(args_dict['log_file'], msg_complete(), 'CONVERT NAMES')

    elif args.cmd == 'normalizeMatrix':
        #Run in sample normalization
        print('Normalizing matrix...')
        run_normalization(args_dict)
        #No os.sys call, no log created when run on own
        #check_process(args_dict['log_file'], msg_complete(), 'NORMALIZE')

    elif args.cmd == 'seRNAseq':
        args_dict['type'] = 'SE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'TRIM')
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        args_dict['input'] = args_dict['alignments']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'ALIGN')
        #Count reads for each alignment file
        msg_count()
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'COUNT')
        #Normalize
        msg_normalize()
        args_dict['input'] = str(args_dict['input']) + str(args_dict['experiment']) + '_counts_table.tsv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'NORMALIZE')
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_metagene(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'QUALITY CONTROL')
        get_summary(args_dict)
        msg_finish()

    elif args.cmd == 'peRNAseq':
        args_dict['type'] = 'PE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'TRIM')
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_peRNAseq(args_dict)
        #Get other formatted files
        args_dict['input'] = args_dict['alignments']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'ALIGN')
        #Count reads for each alignment file
        msg_count()
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'COUNT')
        #Normalize
        msg_normalize()
        args_dict['input'] = str(args_dict['input']) + str(args_dict['experiment']) + '_counts_table.tsv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'NORMALIZE')
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_metagene(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'QUALITY CONTROL')
        get_summary(args_dict)
        msg_finish()

    elif args.cmd == 'riboprof':
        args_dict['type'] = 'SE'
        #Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'TRIM')
        #Align
        msg_align()
        args_dict['input'] = args_dict['trimmed_fastq']
        args_dict = run_seRNAseq(args_dict)
        #Get other formatted files
        args_dict['input'] = args_dict['alignments']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        if args_dict['output_bigwig'] == True:
            create_bigwig(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'ALIGN')
        #Count reads for each alignment file
        msg_count()
        args_dict = count_reads(args_dict)
        #Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'COUNT')
        #Normalize
        msg_normalize()
        args_dict['input'] = str(args_dict['input']) + str(args_dict['experiment']) + '_counts_table.tsv'
        args_dict['gtf'] = args_dict['gtf_type']
        run_normalization(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'NORMALIZE')
        #Run quality control
        msg_quality()
        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)
        args_dict['input'] = args_dict['alignments']
        make_periodicity(args_dict)
        make_metagene(args_dict)
        check_process(args_dict['log_file'], msg_complete(), 'QUALITY CONTROL')
        get_summary(args_dict)
        msg_finish()

    else:
        raise Exception("Invalid function processing function provided.")

"""
DESCRIPTION: Run main
"""
if __name__ == "__main__":

    sys.exit(main() or 0)
