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
import pandas as pd
import datetime

from xpressplot import rpm, r_fpkm, batch_normalize, convert_names #tpm, add back when pip install updated

"""IMPORT INTERNAL DEPENDENCIES"""
from .__init__ import __version__
from .messages import *
from .arguments import get_arguments, get_dependencies
from .trim import run_trim
from .align import run_seRNAseq, run_peRNAseq, create_star_reference
from .buildIndex import index_gtf
from .count import count_reads, collect_counts
from .normalizeMatrix import run_normalization
from .convert import create_bed
from .rrnaProbe import rrnaProbe
from .quality import get_multiqc_summary, get_fastqc
from .metagene import make_metagene
from .geneCoverage import make_coverage
from .readDistribution import make_readDistributions
from .periodicity import make_periodicity
from .complexity import make_complexity
from .parallel import get_cores
from .gtfModify import edit_gtf
from .utils import get_files, get_directories, unzip_files, cleanup
from .test import test_install
from .buildCommand import build_command

"""Main function to call necessary functions for sub-modules

ASSUMPTIONS: Proper arguments are provided where some user renaming of files may be required
"""

"""Wrapper for quantification utilities
"""
def run_quantify(
        args_dict):

    if args_dict['quantification_method'] == 'cufflinks':

        args_dict = count_reads(args_dict)

        # Collect counts into a single table
        args_dict['input'] = args_dict['abundances']
        collect_counts(args_dict)

    elif args_dict['quantification_method'] == 'both':

        args_dict['quantification_method'] = 'htseq'
        args_dict['stranded'] = args_dict['htseq_stranded']
        args_dict = count_reads(args_dict)

        # Collect counts into a single table
        revert_input = args_dict['input']
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)

        args_dict['input'] = revert_input
        args_dict['quantification_method'] = 'cufflinks'
        args_dict['stranded'] = args_dict['cufflinks_stranded']
        args_dict = count_reads(args_dict)

        # Collect counts into a single table
        args_dict['input'] = args_dict['abundances']
        collect_counts(args_dict)

    else:

        args_dict = count_reads(args_dict)

        # Collect counts into a single table
        args_dict['input'] = args_dict['counts']
        collect_counts(args_dict)

    check_process(
        args_dict['log_file'],
        msg_complete(),
        'COUNT') # Check log file for errors and exceptions

    # Normalize
    if args_dict['quantification_method'] == 'htseq':

        msg_normalize()
        args_dict['input'] = str(args_dict['input']) + str(args_dict['experiment']) + '_count_table.tsv'
        run_normalization(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'NORMALIZE') # Check log file for errors and exceptions

    else:
        pass

    return args_dict

# For eventual GUI integration
# from gooey import Gooey
# @Gooey
def main(
    args=None):

    # Read in arguments
    args, args_dict = get_arguments(
        args,
        __version__)

    # Should have already seen check_directory() so should have a trailing '/'
    if 'input' in args_dict \
    and str(args_dict['input']).endswith('/'):
        unzip_files(
            args_dict['input'],
            args_dict['log'])

    # Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'build':
        build_command()

    elif args.cmd == 'test':
        test_install()

    elif args.cmd == 'trim':
        print('Trimming reads...')

        # Trim reads
        run_trim(args_dict)

        print('Generating multiqc summary...')
        get_multiqc_summary(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'TRIM')

    elif args.cmd == 'align':
        print('Aligning reads to reference...')

        # Align
        if args_dict['type'].upper() == 'SE':
            args_dict = run_seRNAseq(args_dict)
        elif args_dict['type'].upper() == 'PE':
            args_dict = run_peRNAseq(args_dict)
        else:
            raise Exception('Invalid type argument provided')

        # Get other formatted files
        args_dict['input'] = args_dict['alignments_coordinates']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)

        print('Generating multiqc summary...')
        get_multiqc_summary(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'ALIGN')

    elif args.cmd == 'count':
        print('Counting alignments...')

        # Count reads for each alignment file
        msg_count()
        args_dict = run_quantify(args_dict)

        print('Generating multiqc summary...')
        get_multiqc_summary(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'COUNT')

    elif args.cmd == 'diffxpress':
        print('Performing differential expression analysis...')

        # Run differential expression analysis via DESeq2
        if str(args_dict['input']).endswith('.txt') or str(args_dict['input']).endswith('.tsv'):
            output_file = str(args_dict['input'])[:-4] + '_diffx.tsv'
        else:
            raise Exception('Unrecognized input_file delimiter type. Files must be tab-delimited')

        if str(args_dict['sample']).endswith('.txt') or str(args_dict['sample']).endswith('.tsv'):
            pass
        else:
            raise Exception('Unrecognized sample_file delimiter type. Files must be tab-delimited')

        if str(args_dict['design']).startswith('~'):
            raise Exception('Tilde should not be included in design formula, script will automatically add this syntax.')

        if 'shrink' in args_dict \
        and args_dict['shrink'] == True:
            shrinkage = 'TRUE'
        else:
            shrinkage = 'FALSE'

        # Run deseq2 in R
        os.system('Rscript' \
            + ' ' + str(args_dict['path']) + 'Rdiffxpress.r' \
            + ' ' + str(args_dict['input']) \
            + ' ' + str(args_dict['sample']) \
            + ' ' + str(output_file) \
            + ' ' + str(args_dict['design'])
            + ' ' + str(shrinkage)
            + str(args_dict['log']))

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'DIFFERENTIAL EXPRESSION')

    elif args.cmd == 'metagene':
        print('Performing metagene analysis on transcriptome-mapped files...')

        # Get list of bam files from user input
        files = get_files(
            args_dict['input'],
            [str(args_dict['bam_suffix'])])
        if len(files) == 0:
            raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

        # Perform metagene analysis
        args_dict['gene_name'] = None
        index_gtf(args_dict)
        make_metagene(args_dict, files)

        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.idx')

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'METAGENE')

    elif args.cmd == 'geneCoverage':
        print('Performing gene coverage analysis on transcriptome-mapped files...')

        # Perform metagene analysis
        # Get list of bam files from user input
        files = get_files(
            args_dict['input'],
            [str(args_dict['bam_suffix'])])
        if len(files) == 0:
            raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

        success = index_gtf(args_dict, gene_name=args_dict['gene_name'])

        if success != -1:
            make_coverage(args_dict, files)
            os.system(
                'rm'
                + ' ' + args_dict['output'] + '*.fts')
        else:
            print('Could not find ' + str(args_dict['gene_name']) + ' in reference. Please try running the geneCoverage module with another known gene for your organism.')

        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.idx')
        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.fts')

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'GENE COVERAGE')

    elif args.cmd == 'readDistribution':
        print('Performing read distribution analysis on fastq files...')

        # Generate read distribution summaries
        make_readDistributions(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'READ DISTRIBUTION')

    elif args.cmd == 'p_sites':
        print('Performing P-site analysis on transcriptome-mapped files...')

        # Generate P-site summaries
        make_periodicity(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'PERIODICITY')

    elif args.cmd == 'complexity':
        print('Performing library complexity analysis...')

        # Generate library complexity summaries
        make_complexity(args_dict)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'COMPLEXITY')

    elif args.cmd == 'curateReference':
        print('Curating reference')

        # Create STAR reference
        args_dict['threads'], args_dict['workers'] = get_cores(
            args_dict,
            mod_workers = True)
        create_star_reference(
            args_dict['output'],
            args_dict['fasta'],
            args_dict['gtf'],
            args_dict['log'],
            threads = args_dict['threads'],
            sjdbOverhang = args_dict['sjdbOverhang'],
            genome_size = args_dict['genome_size'])

        # Truncate transcript reference
        if (args_dict['longest_transcript'] == True) or (args_dict['protein_coding'] == True) or (args_dict['truncate'] == True):
            edit_gtf(
                args_dict['gtf'],
                longest_transcript = args_dict['longest_transcript'],
                protein_coding = args_dict['protein_coding'],
                truncate_reference = args_dict['truncate'],
                _5prime = args_dict['truncate_5prime'], # If no 5' truncation desired, set to 0
                _3prime = args_dict['truncate_3prime'], # If no 3' truncation desired, set to 0
                output = True,
                threads = args_dict['threads'])

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'CURATE REFERENCE')

    elif args.cmd == 'makeReference':
        print('Creating reference files...')

        # Generate reference
        args_dict['threads'], args_dict['workers'] = get_cores(
            args_dict,
            mod_workers = True)
        create_star_reference(
            args_dict['output'],
            args_dict['fasta'],
            args_dict['gtf'],
            args_dict['log'],
            threads = args_dict['threads'],
            sjdbOverhang = args_dict['sjdbOverhang'],
            genome_size = args_dict['genome_size'])

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'MAKE REFERENCE')

    elif args.cmd == 'modifyGTF':
        print('Formatting reference file...')

        # Truncate transcript reference
        args_dict['threads'], args_dict['workers'] = get_cores(
            args_dict,
            mod_workers=True)
        edit_gtf(
            args_dict['gtf'],
            longest_transcript = args_dict['longest_transcript'],
            protein_coding = args_dict['protein_coding'],
            truncate_reference = args_dict['truncate'],
            _5prime = args_dict['truncate_5prime'], # If no 5' truncation desired, set to 0
            _3prime = args_dict['truncate_3prime'], # If no 3' truncation desired, set to 0
            output = True,
            threads = args_dict['threads'])

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        # check_process(args_dict['log_file'], msg_complete(), 'TRUNCATE')

    elif args.cmd == 'rrnaProbe':

        # Get files to probe
        probe_list = get_directories(
            args_dict['input'],
            [''],
            omit=['.html','.zip'])

        # Run rrna_prober, output to outputDir
        print('Probing for most over-represented read sequences...')
        probe_out = rrnaProbe(
            probe_list,
            args_dict['min_overlap']) # Use inputDir to get FASTQC files and output to outputDir/analysis

        # Output summary
        with open(args_dict['output'] + 'rrnaProbe_output.txt', "w") as text_file:
            print(probe_out, file = text_file)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'RRNA PROBE')

    elif args.cmd == 'convertNames':
        print('Converting row names...')

        # Convert row names in dataframe
        if str(args_dict['input']).endswith('.csv'):
            delim = ','
            suf = '.csv'
        else:
            delim = '\t'
            suf = '.tsv'

        data = pd.read_csv(str(args_dict['input']), sep=delim, index_col=0)
        data = convert_names(
            data,
            args_dict['gtf'],
            orig_name_label = args_dict['orig_name_label'],
            orig_name_location = args_dict['orig_name_location'],
            new_name_label = args_dict['new_name_label'],
            new_name_location = args_dict['new_name_location'],
            refill = args_dict['refill'],
            sep = '\t')

        data.to_csv(
            str(args_dict['input'])[:-4] + '_renamed' + str(suf),
            sep = delim)

        # Check log file for errors and exceptions
        get_dependencies(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'CONVERT NAMES')

    elif args.cmd == 'normalizeMatrix':
        print('Normalizing matrix...')

        #Run in sample normalization
        run_normalization(args_dict)
        get_dependencies(args_dict)

        # Check log file for errors and exceptions
        #No os.sys call, no log created when run on own
        #check_process(args_dict['log_file'], msg_complete(), 'NORMALIZE')

    elif args.cmd == 'seRNAseq':
        args_dict['type'] = 'SE'

        # Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'TRIM') # Check log file for errors and exceptions

        # Align
        args_dict['input'] = args_dict['trimmed_fastq']
        msg_fastqc()
        get_fastqc(args_dict) # Run FastQC on trimmed reads
        msg_align()
        args_dict = run_seRNAseq(args_dict)

        # Get other formatted files
        args_dict['input'] = args_dict['alignments_coordinates']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'ALIGN') # Check log file for errors and exceptions

        # Count reads for each alignment file
        msg_count()
        args_dict = run_quantify(args_dict)

        # Run quality control
        msg_quality()

        # Get multiqc report and print close message
        get_multiqc_summary(args_dict)

        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)

        args_dict['input'] = args_dict['alignments_coordinates']
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        make_complexity(args_dict)

        args_dict['input'] = args_dict['alignments_transcriptome']
        args_dict['bam_suffix'] = 'toTranscriptome.out.bam'
        args_dict['gene_name'] = 'GAPDH'
        args_dict['samples'] = None

        # Get list of bam files from user input
        files = get_files(
            args_dict['input'],
            [str(args_dict['bam_suffix'])])
        if len(files) == 0:
            raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

        success = index_gtf(args_dict, gene_name=args_dict['gene_name'])

        if success != -1:
            make_coverage(args_dict, files)
            os.system(
                'rm'
                + ' ' + args_dict['output'] + '*.fts')
        else:
            print('Could not find ' + str(args_dict['gene_name']) + ' in reference. Please try running the geneCoverage module with another known housekeeping gene for your organism.')

        args_dict['gene_name'] = None
        index_gtf(args_dict)
        make_metagene(args_dict, files)
        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.idx')

        if 'small_output' in args_dict and args_dict['small_output'] == True:
            cleanup(args_dict)

        check_process(
            args_dict['log_file'],
            msg_complete(),
            'QUALITY CONTROL') # Check log file for errors and exceptions

        get_dependencies(args_dict)
        msg_finish()

    elif args.cmd == 'peRNAseq':
        args_dict['type'] = 'PE'

        # Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'TRIM') # Check log file for errors and exceptions

        # Align
        args_dict['input'] = args_dict['trimmed_fastq']
        msg_fastqc()
        get_fastqc(args_dict) # Run FastQC on trimmed reads
        msg_align()
        args_dict = run_peRNAseq(args_dict)

        # Get other formatted files
        args_dict['input'] = args_dict['alignments_coordinates']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'ALIGN') # Check log file for errors and exceptions

        # Count reads for each alignment file
        msg_count()
        args_dict = run_quantify(args_dict)

        # Run quality control
        msg_quality()
        # Get multiqc report and print close message
        get_multiqc_summary(args_dict)

        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)

        args_dict['input'] = args_dict['alignments_coordinates']
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        make_complexity(args_dict)

        args_dict['input'] = args_dict['alignments_transcriptome']
        args_dict['bam_suffix'] = 'toTranscriptome.out.bam'
        args_dict['gene_name'] = 'GAPDH'
        args_dict['samples'] = None

        # Get list of bam files from user input
        files = get_files(
            args_dict['input'],
            [str(args_dict['bam_suffix'])])
        if len(files) == 0:
            raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

        success = index_gtf(args_dict, gene_name=args_dict['gene_name'])

        if success != -1:
            make_coverage(args_dict, files)
            os.system(
                'rm'
                + ' ' + args_dict['output'] + '*.fts')
        else:
            print('Could not find ' + str(args_dict['gene_name']) + ' in reference. Please try running the geneCoverage module with another known housekeeping gene for your organism.')

        args_dict['gene_name'] = None
        index_gtf(args_dict)
        make_metagene(args_dict, files)

        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.idx')

        if 'small_output' in args_dict and args_dict['small_output'] == True:
            cleanup(args_dict)

        check_process(
            args_dict['log_file'],
            msg_complete(),
            'QUALITY CONTROL') # Check log file for errors and exceptions

        get_dependencies(args_dict)
        msg_finish()

    elif args.cmd == 'riboseq':
        args_dict['type'] = 'SE'

        # Trim
        msg_trim()
        args_dict = run_trim(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'TRIM') # Check log file for errors and exceptions

        # Align
        args_dict['input'] = args_dict['trimmed_fastq']
        msg_fastqc()
        get_fastqc(args_dict) # Run FastQC on trimmed reads
        msg_align()
        args_dict = run_seRNAseq(args_dict)

        # Get other formatted files
        args_dict['input'] = args_dict['alignments_coordinates']
        if args_dict['output_bed'] == True:
            create_bed(args_dict)
        check_process(
            args_dict['log_file'],
            msg_complete(),
            'ALIGN') # Check log file for errors and exceptions

        # Count reads for each alignment file
        msg_count()
        args_dict = run_quantify(args_dict)

        # Run quality control
        msg_quality()

        # Get multiqc report and print close message
        get_multiqc_summary(args_dict)

        args_dict['input'] = args_dict['trimmed_fastq']
        make_readDistributions(args_dict)

        args_dict['input'] = args_dict['alignments_coordinates']
        args_dict['gtf'] = str(args_dict['reference']) + 'transcripts.gtf'
        make_complexity(args_dict)

        args_dict['input'] = args_dict['alignments_transcriptome']
        args_dict['bam_suffix'] = 'toTranscriptome.out.bam'
        if 'gene_name' not in args_dict \
        or args_dict['gene_name'] == None:
            args_dict['gene_name'] = 'GAPDH'
        args_dict['samples'] = None

        # Get list of bam files from user input
        files = get_files(
            args_dict['input'],
            [str(args_dict['bam_suffix'])])
        if len(files) == 0:
            raise Exception('No files with suffix ' + str(args_dict['bam_suffix']) + ' found in the directory ' +  str(args_dict['input']))

        success = index_gtf(args_dict, gene_name=args_dict['gene_name'])

        if success != -1:

            make_coverage(args_dict, files)
            os.system(
                'rm'
                + ' ' + args_dict['output'] + '*.fts')
        else:
            print('Could not find ' + str(args_dict['gene_name']) + ' in reference. Please try running the geneCoverage module with another known housekeeping gene for your organism.')

        args_dict['gene_name'] = None
        index_gtf(args_dict)
        make_metagene(args_dict, files)
        make_periodicity(args_dict)

        os.system(
            'rm'
            + ' ' + args_dict['output'] + '*.idx')

        if 'small_output' in args_dict and args_dict['small_output'] == True:
            cleanup(args_dict)

        check_process(
            args_dict['log_file'],
            msg_complete(),
            'QUALITY CONTROL') # Check log file for errors and exceptions

        get_dependencies(args_dict)
        msg_finish()

    else:
        raise Exception('Invalid function processing function provided.')

"""Run main"""
if __name__ == '__main__':

    sys.exit(main() or 0)
