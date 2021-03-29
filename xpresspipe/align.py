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
import re
import gc
import csv
import pandas as pd

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, add_directory, get_fasta, check_directories
#from .utils import clean_vcf
from .gtfModify import edit_gtf
from .parallel import parallelize, parallelize_pe

gtf_type_column = 2
gtf_annotation_column = 8
search_type = 'transcript'
parse_type = 'rrna'

def generate_bed(
        gtf_file):
    """Generate a BED file of rRNA sequences for depletion from genome-aligned
    BAM file.
    """

    # Set up BED generation command
    gtf = pd.read_csv(
        str(gtf_file),
        sep = '\t',
        header = None,
        comment = '#',
        low_memory = False)

    # Remove records that map to rRNA
    gtf_rrna = gtf[gtf[gtf_annotation_column].str.contains(
        parse_type, flags=re.IGNORECASE
        )]
    gtf_rrna = gtf_rrna.loc[gtf_rrna[gtf_type_column] == search_type]
    gtf_rrna = gtf_rrna[[0,3,4]]

    bed_file = gtf_file[:-4] + '_rrna.bed'

    gtf_rrna.to_csv(
        str(bed_file),
        sep = '\t',
        header = None,
        index = False,
        quoting = csv.QUOTE_NONE)

    gtf = None # Garbage management
    gtf_rrna = None
    gc.collect()

    return bed_file

def curate_depletion_gtf(args_dict):

    gtf = str(args_dict['reference']) + 'transcripts.gtf'

    gtf = pd.read_csv(
        str(gtf),
        sep = '\t',
        header = None,
        comment = '#',
        low_memory = False)

    # Remove records that map to rRNA
    gtf_depl = gtf[~gtf[gtf_annotation_column].str.contains(
        parse_type, flags=re.IGNORECASE
        )]

    gtf_depl.to_csv(
        str(args_dict['reference']) + 'transcripts_depl.gtf',
        sep = '\t',
        header = None,
        index = False,
        quoting = csv.QUOTE_NONE)

    gtf = None # Garbage management
    gtf_depl = None
    gc.collect()

    return 'transcripts_depl.gtf'

"""Create STAR reference genome"""
def create_star_reference(
    output_directory,
    fasta_directory,
    gtf,
    log,
    threads=1,
    sjdbOverhang=100,
    genome_size=14):

    # Create output directory
    output_directory = check_directories(output_directory)
    fasta_directory = check_directories(fasta_directory)

    os.system('mkdir'
        + ' ' + str(output_directory) + 'genome'
        + str(log))

    fasta_list = get_fasta(fasta_directory)

    # Create reference
    os.system(
        'STAR'
        + ' --runMode genomeGenerate'
        + ' --genomeDir ' + str(output_directory) + 'genome'
        + ' --genomeFastaFiles ' + str(fasta_list)
        + ' --sjdbOverhang ' + str(sjdbOverhang)
        + ' --genomeSAindexNbases ' + str(genome_size)
        + ' --sjdbGTFfile ' + str(gtf)
        + ' --runThreadN ' + str(threads)
        + str(log))

"""Build intermediate STAR alignment reference using splice junction annotations from first pass"""
def build_star_splice_junction_intermediate(
    output,
    args_dict):

    os.system(
        'mkdir'
        + ' '
        + str(args_dict['intermediate_references']) + str(output) # Make directory for currently processed file
        + str(args_dict['log'])) # Record log output (must go last in command)

    fasta_list = get_fasta(args_dict['reference'])

    os.system(
        'STAR'
        + ' --runMode genomeGenerate'
        + ' --runThreadN ' + str(args_dict['threads'])
        + ' --sjdbFileChrStartEnd ' + str(args_dict['alignments_coordinates']) + str(output) + '_SJ.out.tab' # Splice junction annotation file output in first pass alignment
        + ' --genomeFastaFiles ' + str(args_dict['fasta_list']) # Input chromosomal fasta files for reference building
        + ' --genomeSAindexNbases ' + str(args_dict['genome_size'])
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output) # Location for output revised reference
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) # Read overhand amount to allow for splice mapping (should be same used in curation of reference)
        + str(args_dict['log'])) # Record log output (must go last in command)

"""Run first pass STAR alignment to map splice junctions"""
def first_pass_star(
    file,
    output,
    args_dict):

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads']) # Argument to specify number of threads to use for processing
        + ' --genomeDir ' + str(args_dict['reference']) + 'genome' # Argument for specifying STAR reference directory
        + ' --readFilesIn ' + str(file) # Argument to dictate directory where pre-processed read files are located
        + ' --outFileNamePrefix ' + str(args_dict['alignments_coordinates']) + str(output) + '_' # Argument to dictate output directory
        + ' --outFilterMismatchNoverLmax ' + str(args_dict['mismatchRatio']) # Mismatch ratio to mapped read length
        + ' --seedSearchStartLmax ' + str(args_dict['seedSearchStartLmax']) #
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) # Read overhand amount to allow for splice mapping (should be same used in curation of reference)
        + ' --outFilterMultimapScoreRange 1'
        + ' --outFilterMultimapNmax 20'
        + ' --outFilterMismatchNmax 10'
        + ' --alignIntronMax 500000'
        + ' --alignMatesGapMax 1000000'
        + ' --sjdbScore 2'
        + ' --alignSJDBoverhangMin 1'
        + ' --genomeLoad NoSharedMemory'
        + ' --readFilesCommand cat'
        + ' --outFilterMatchNminOverLread 0.33'
        + ' --outFilterScoreMinOverLread 0.33'
        + ' --outSAMstrandField intronMotif'
        + ' --outSAMtype None'
        + ' --outSAMmode None'
        + str(args_dict['log'])) # Record log output (must go last in command)

"""Run second pass STAR alignment to map reads splice-aware"""
def second_pass_star(
    file,
    output,
    args_dict):

    if 'vcf' in args_dict and args_dict['vcf'] != None:
        #clean_vcf(args_dict)
        #args_dict['vcf'] = str(args_dict['reference']) + 'xpressyourself.vcf'
        attr_cols = ' --outSAMattributes NH HI NM MD AS XS vG vA'
        vcf_line = ' --varVCFfile ' + str(args_dict['vcf'])
    else:
        attr_cols = ' --outSAMattributes NH HI NM MD AS XS'
        vcf_line = ''

    if args_dict['remove_rrna'] == True:
        gtf_name = curate_depletion_gtf(args_dict)
    else:
        gtf_name = 'transcripts.gtf'

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads'])
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output)
        + ' --readFilesIn ' + str(file)
        + ' --sjdbGTFfile ' + str(args_dict['reference']) + str(gtf_name)
        + ' --outFileNamePrefix ' + str(args_dict['alignments_coordinates']) + str(output) + '_'
        + ' --outFilterMismatchNoverLmax ' + str(args_dict['mismatchRatio']) # Mismatch ratio to mapped read length
        + ' --seedSearchStartLmax ' + str(args_dict['seedSearchStartLmax'])
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang'])
        + ' --outFilterMultimapScoreRange 1'
        + ' --outFilterMultimapNmax 20'
        + ' --outFilterMismatchNmax 10'
        + ' --alignIntronMax 500000'
        + ' --alignMatesGapMax 1000000'
        + ' --sjdbScore 2'
        + ' --alignSJDBoverhangMin 1'
        + ' --genomeLoad NoSharedMemory'
        + ' --limitBAMsortRAM 0'
        + ' --readFilesCommand cat'
        + ' --outFilterMatchNminOverLread 0.33'
        + ' --outFilterScoreMinOverLread 0.33'
        + ' --outSAMstrandField intronMotif'
        + str(attr_cols)
        + str(vcf_line)
        + ' --outSAMunmapped Within'
        + ' --outSAMtype BAM Unsorted' # Allow for multithreading STAR run without file overload
        + ' --quantMode TranscriptomeSAM'
        + ' --outSAMheaderHD @HD VN:1.4'
        + str(args_dict['log']))

"""Run one-pass STAR with GTF reference annotation to guide mapping to splice regions"""
def guided_star(
    file,
    output,
    args_dict):

    if 'vcf' in args_dict and args_dict['vcf'] != None:
        #clean_vcf(args_dict)
        #args_dict['vcf'] = str(args_dict['reference']) + 'xpressyourself.vcf'
        attr_cols = ' --outSAMattributes NH HI NM MD AS XS vG vA'
        vcf_line = ' --varVCFfile ' + str(args_dict['vcf'])
    else:
        attr_cols = ' --outSAMattributes NH HI NM MD AS XS'
        vcf_line = ''

    if args_dict['remove_rrna'] == True:
        gtf_name = curate_depletion_gtf(args_dict)
    else:
        gtf_name = 'transcripts.gtf'

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads']) # Argument to specify number of threads to use for processing
        + ' --genomeDir ' + str(args_dict['reference']) + 'genome' # Argument for specifying STAR reference directory
        + ' --readFilesIn ' + str(file) # Argument to dictate directory where pre-processed read files are located
        + ' --sjdbGTFfile ' + str(args_dict['reference']) + str(gtf_name)
        + ' --outFileNamePrefix ' + str(args_dict['alignments_coordinates']) + str(output) + '_'
        + ' --outFilterMismatchNoverLmax ' + str(args_dict['mismatchRatio']) # Mismatch ratio to mapped read length
        + ' --seedSearchStartLmax ' + str(args_dict['seedSearchStartLmax'])
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang'])
        + ' --outFilterMultimapScoreRange 1'
        + ' --outFilterMultimapNmax 20'
        + ' --outFilterMismatchNmax 10'
        + ' --alignIntronMax 500000'
        + ' --alignMatesGapMax 1000000'
        + ' --sjdbScore 2'
        + ' --alignSJDBoverhangMin 1'
        + ' --genomeLoad NoSharedMemory'
        + ' --limitBAMsortRAM 0'
        + ' --readFilesCommand cat'
        + ' --outFilterMatchNminOverLread 0.33'
        + ' --outFilterScoreMinOverLread 0.33'
        + ' --outSAMstrandField intronMotif'
        + str(attr_cols)
        + str(vcf_line)
        + ' --outSAMunmapped Within'
        + ' --outSAMtype BAM Unsorted' # Allow for multithreading STAR run without file overload
        + ' --quantMode TranscriptomeSAM'
        + ' --outSAMheaderHD @HD VN:1.4'
        + str(args_dict['log']))

"""Remove intermediate reference files for the file being processed in the instance"""
def remove_intermediate_reference(
    output,
    args_dict):

    # Clear the current file's splice junction intermediate reference
    os.system(
        'rm -r'
        + ' ' + str(args_dict['intermediate_references']) + str(output)
        + str(args_dict['log']))

"""Sort reads per file by chromosome position and keep only unique mappers"""
def alignment_process(
    output,
    args_dict,
    paired=False):

    # Fixmates for paired-end
    if paired == True:

        # Fixmate requires name sorted files
        os.system(
        'samtools sort'
        + ' -n'
        + ' --threads ' + str(args_dict['threads'])
        + ' -o ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.namesort.bam'
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.out.bam')

        # Run fixmate
        os.system(
        'samtools fixmate'
        + ' -m'
        + ' --threads ' + str(args_dict['threads'])
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.namesort.bam'
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_fixed.namesort.bam')

        # Remove rRNA
        if args_dict['remove_rrna'] == True:
            os.system(
                'bedtools intersect -abam'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_fixed.namesort.bam'
                + ' -b ' + args_dict['bed_file']
                + ' -v > ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_fixed.rrna_depl.bam')
            os.system(
                'mv'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_fixed.rrna_depl.bam'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_fixed.namesort.bam')

        # Convert back to coordinate sorted because markdup doesn't accept name sorted files
        os.system(
        'samtools sort'
        + ' --threads ' + str(args_dict['threads'])
        + ' -o ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.sort.bam'
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_fixed.namesort.bam')

    else:
        # Remove rRNA
        if args_dict['remove_rrna'] == True:
            os.system(
                'bedtools intersect -abam'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_Aligned.out.bam'
                + ' -b ' + args_dict['bed_file']
                + ' -v > ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_Aligned.rrna_depl.bam')
            os.system(
                'mv'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_Aligned.rrna_depl.bam'
                + ' ' + str(args_dict['alignments_coordinates'])
                + str(output) + '_Aligned.out.bam')

        # Sort SAM file
        os.system(
            'samtools sort'
            + ' --threads ' + str(args_dict['threads'])
            + ' -o ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.sort.bam'
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.out.bam'
            + str(args_dict['log']))

    # Only take unique mappers (q = 255)
    if args_dict['no_multimappers'] == True:
        os.system(
            'samtools view'
            + ' -h' # Keep SAM header in output
            + ' -q 255' # Keep unique mappers
            + ' --threads ' + str(args_dict['threads'])
            + ' -o ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.unique.bam'
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.sort.bam'
            + str(args_dict['log']))
        os.system(
            'mv '
            + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.unique.bam'
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_Aligned.sort.bam')
        file_suffix = '_Aligned.sort.bam'
    else:
        file_suffix = '_Aligned.sort.bam'

    # Index BAM file
    os.system(
        'samtools index'
        + ' -@ ' + str(args_dict['threads'])
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + str(file_suffix)
        + str(args_dict['log']))

    if ('umi_location' in args_dict
    and (args_dict['umi_location'] != None
    and str(args_dict['umi_location']).lower() != 'none')
    ) or ('umi_seq' in args_dict and args_dict['umi_seq'] == True):
        os.system(
            'umi_tools dedup'
            + ' -I ' + str(args_dict['alignments_coordinates']) + str(output) + str(file_suffix) # Input BAM
            + ' -S ' + str(args_dict['alignments_coordinates']) + str(output) + '_UMIremoved.bam' # Output BAM
            + str(args_dict['log']))
        os.system(
            'samtools index'
            + ' -@ ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_UMIremoved.bam'
            + str(args_dict['log']))

        # Use sorted BAM file to find any duplicate reads
        os.system(
            'samtools markdup'
            + ' --threads ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_UMIremoved.bam' # Input BAM
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_UMImarked.bam' # Output BAM
            + ' -s' # Print some basic stats
            + str(args_dict['log']))
        os.system(
            'samtools index'
            + ' -@ ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_UMImarked.bam'
            + str(args_dict['log']))
    else:
        # Use sorted BAM file to find any duplicate reads
        os.system(
            'samtools markdup'
            + ' --threads ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + str(file_suffix) # Input BAM
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_dedupMarked.bam' # Output BAM
            + ' -s' # Print some basic stats
            + str(args_dict['log']))
        os.system(
            'samtools index'
            + ' -@ ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_dedupMarked.bam'
            + str(args_dict['log']))

        # Create sorted BAM file with duplicates removed
        os.system(
            'samtools markdup'
            + ' --threads ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + str(file_suffix) # Input BAM
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_dedupRemoved.bam' # Output BAM
            + ' -s' # Print some basic stats to STDOUT
            + ' -r' # Remove duplicate reads
            + str(args_dict['log']))
        os.system(
            'samtools index'
            + ' -@ ' + str(args_dict['threads'])
            + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '_dedupRemoved.bam'
            + str(args_dict['log']))

def store_alignments_transcriptome(
    output,
    args_dict):

    os.system(
        'mv'
        + ' ' + str(args_dict['alignments_coordinates']) + str(output) + '*Aligned.toTranscriptome.out.bam'
        + ' ' + str(args_dict['output']) + 'alignments_transcriptome' )

"""Remove all intermediate alignment files and references after alignment is complete"""
def remove_intermediates(
    args_dict):

    os.system(
        "find"
        + " " + str(args_dict['alignments_coordinates'])
        + " -maxdepth 1 -type f" # Only keep files matching pattern
        + " ! -name '*_Aligned.sort.bam'"
        + " ! -name '*_Aligned.sort.bam.bai'"
        + " ! -name '*_fixed.sort.bam'"
        + " ! -name '*_fixed.sort.bam.bai'"
        + " ! -name '*_dedupMarked.bam'"
        + " ! -name '*_dedupMarked.bam.bai'"
        + " ! -name '*_dedupRemoved.bam'"
        + " ! -name '*_dedupRemoved.bam.bai'"
        + " ! -name '*_UMIremoved.bam'"
        + " ! -name '*_UMIremoved.bam.bai'"
        + " ! -name '*_UMImarked.bam'"
        + " ! -name '*_UMImarked.bam.bai'"
        + " ! -name '*_Log.final.out'"
        + " -delete"
        + str(args_dict['log']))

def clean_reference_directory(
    args_dict):

    os.system(
        'rm -r'
        + ' ' + str(args_dict['alignments_coordinates']) + '*_STARgenome'
        + str(args_dict['log']))

    if 'two-pass' in args_dict and args_dict['two-pass'] == True:
        os.system(
            'rm -r'
            + ' ' + str(args_dict['intermediate_references'])
            + str(args_dict['log']))

def align(
    args_dict,
    output,
    file,
    paired=False):

    if args_dict['threads'] > 30:
        thread_count = args_dict['threads']
        args_dict['threads'] = 30
    else:
        thread_count = args_dict['threads']

    #if 'mask' in args_dict and args_dict['mask'] == True:
    #    file = masking_star(
    #        file,
    #        output,
    #        args_dict)

    if 'two-pass' in args_dict and args_dict['two-pass'] == True:

        # STAR first pass
        first_pass_star(
            file,
            output,
            args_dict)

        # STAR intermediate reference building
        build_star_splice_junction_intermediate(
            output,
            args_dict)

        # STAR second pass
        second_pass_star(
            file,
            output,
            args_dict)

        # Remove intermediate reference for the file
        remove_intermediate_reference(
            output,
            args_dict)

    else:
        # One-pass STAR with GTF guiding splice mapping
        guided_star(
            file,
            output,
            args_dict)
        store_alignments_transcriptome(
            output,
            args_dict)

    # Create BAM file with only unique hits, mark duplicates, index
    alignment_process(
        output,
        args_dict,
        paired=paired)

    # Clean up the output
    remove_intermediates(
        args_dict)

    #if 'vcf' in args_dict and args_dict['vcf'] != None:
    #    os.system(
    #        'rm'
    #        + ' ' + str(args_dict['reference']) + 'xpressyourself.vcf')

    args_dict['threads'] = thread_count

"""Single-end RNA-seq pipeline"""
def se_align(
    args):

    file, args_dict = args[0], args[1]

    output = str(file).rsplit('.',1)[0].replace('trimmed_','')
    file = str(args_dict['input']) + str(file)
    align(args_dict, output, file)

"""Paired-end RNA-seq pipeline"""
def pe_align(
    args):

    file1, file2, args_dict = args[0], args[1], args[2]

    # STAR first pass
    file = str(args_dict['input']) + str(file1) + ' ' + str(args_dict['input']) + str(file2)
    output = str(file1).rsplit('.',1)[0].replace('trimmed_','').replace('read1','').replace('read2','').replace('r1','').replace('r2','')

    align(args_dict, output, file, paired=True)

"""Manage single-end RNA-seq pipeline"""
def run_seRNAseq(
    args_dict):

    # Add output directories
    args_dict = add_directory(
        args_dict,
        'output',
        'alignments_coordinates')
    args_dict = add_directory(
        args_dict,
        'output',
        'alignments_transcriptome')

    if args_dict['remove_rrna'] == True:
        gtf_file = str(args_dict['reference']) + 'transcripts.gtf'
        args_dict['bed_file'] = generate_bed(
            gtf_file=gtf_file)

    if 'two-pass' in args_dict and args_dict['two-pass'] == True:
        args_dict = add_directory(
            args_dict,
            'alignments_coordinates',
            'intermediate_references')

    # Get list of files to align based on acceptable file types
    files = get_files(
        args_dict['input'],
        ['.fastq','.fq','.txt'])

    # Align single-end RNAseq reads
    parallelize(
        se_align,
        files,
        args_dict)
    clean_reference_directory(args_dict)

    return args_dict


"""Manage paired-end RNA-seq pipeline"""
def run_peRNAseq(
    args_dict):

    # Add output directories
    args_dict = add_directory(
        args_dict,
        'output',
        'alignments_coordinates')
    args_dict = add_directory(
        args_dict,
        'output',
        'alignments_transcriptome')

    if args_dict['remove_rrna'] == True:
        gtf_file = str(args_dict['reference']) + 'transcripts.gtf'
        args_dict['bed_file'] = generate_bed(
            gtf_file=gtf_file)

    if 'two-pass' in args_dict and args_dict['two-pass'] == True:
        args_dict = add_directory(
            args_dict,
            'alignments_coordinates',
            'intermediate_references')

    # Get list of files to align based on acceptable file types
    files = get_files(
        args_dict['input'],
        ['.fastq','.fq','.txt'])

    if len(files) % 2 != 0:
        raise Exception('An uneven number of paired-end files were specified in the input directory')
    else:
        # Align paired-end RNAseq reads
        parallelize_pe(
            pe_align,
            files,
            args_dict)
        clean_reference_directory(args_dict)

    return args_dict
