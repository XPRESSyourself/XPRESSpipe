"""XPRESSpipe
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

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, add_directory, get_fasta, check_directories
from .parallel import parallelize, parallelize_pe

"""Create STAR reference genome"""
def create_star_reference(
        output_directory,
        fasta_directory,
        gtf,
        log,
        threads=1,
        sjdbOverhang=100):

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
        + ' --sjdbGTFfile ' + str(gtf)
        + ' --runThreadN ' + str(threads)
        + str(log))

"""Create STAR mask reference"""
def create_mask_reference(
        output_directory,
        fasta_directory,
        log,
        threads=1):

    # Create output directory
    output_directory = check_directories(output_directory)
    fasta_directory = check_directories(fasta_directory)

    os.system('mkdir'
        + ' ' + str(output_directory) + 'mask'
        + str(log))

    fasta_list = get_fasta(fasta_directory)

    # Create reference
    os.system(
        'STAR'
        + ' --runMode genomeGenerate'
        + ' --genomeDir ' + str(output_directory) + 'mask'
        + ' --genomeFastaFiles ' + str(fasta_list)
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

    os.system(
        'STAR'
        + ' --runMode genomeGenerate'
        + ' --runThreadN ' + str(args_dict['threads'])
        + ' --sjdbFileChrStartEnd ' + str(args_dict['alignments']) + str(output) + '_SJ.out.tab' # Splice junction annotation file output in first pass alignment
        + ' --genomeFastaFiles ' + str(args_dict['fasta_list']) # Input chromosomal fasta files for reference building
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output) # Location for output revised reference
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) # Read overhand amount to allow for splice mapping (should be same used in curation of reference)
        + str(args_dict['log'])) # Record log output (must go last in command)

"""Run a masking alignment step with STAR"""
def masking_star(
    file,
    output,
    args_dict):

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads']) # Argument to specify number of threads to use for processing
        + ' --genomeDir ' + str(args_dict['reference']) + 'mask' # Argument for specifying STAR reference directory
        + ' --readFilesIn ' + str(file) # Argument to dictate directory where pre-processed read files are located
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_' # Argument to dictate output directory
        + ' --outReadsUnmapped Fastx'
        + ' --outFilterMultimapScoreRange 1'
        + ' --outFilterMultimapNmax 20'
        + ' --outFilterMismatchNmax 10'
        + ' --alignIntronMax 500000'
        + ' --alignMatesGapMax 1000000'
        + ' --genomeLoad NoSharedMemory'
        + ' --readFilesCommand cat'
        + ' --outFilterMatchNminOverLread 0.33'
        + ' --outFilterScoreMinOverLread 0.33'
        + str(args_dict['log']))

    if ' ' in file:
        return str(args_dict['alignments']) + str(output) + '_Unmapped.out.mate1 ' + str(args_dict['alignments']) + str(output) + '_Unmapped.out.mate2'
    else:
        return str(args_dict['alignments']) + str(output) + '_Unmapped.out.mate1'

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
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_' # Argument to dictate output directory
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

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads'])
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output)
        + ' --readFilesIn ' + str(file)
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_'
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
        + ' --outSAMattributes NH HI NM MD AS XS'
        + ' --outSAMunmapped Within'
        + ' --outSAMtype BAM Unsorted' # Allow for multithreading STAR run without file overload
        + ' --outSAMheaderHD @HD VN:1.4'
        + str(args_dict['log']))

"""Run one-pass STAR with GTF reference annotation to guide mapping to splice regions"""
def guided_star(
    file,
    output,
    args_dict):

    os.system(
        'STAR'
        + ' --runThreadN ' + str(args_dict['threads']) # Argument to specify number of threads to use for processing
        + ' --genomeDir ' + str(args_dict['reference']) + 'genome' # Argument for specifying STAR reference directory
        + ' --readFilesIn ' + str(file) # Argument to dictate directory where pre-processed read files are located
        + ' --sjdbGTFfile ' + str(args_dict['reference']) + 'transcripts.gtf'
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_'
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
        + ' --outSAMattributes NH HI NM MD AS XS'
        + ' --outSAMunmapped Within'
        + ' --outSAMtype BAM Unsorted' # Allow for multithreading STAR run without file overload
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
        args_dict):

    # Only take unique mappers (q = 255)
    os.system(
        'samtools view'
        + ' -h' # Keep SAM header in output
        + ' -q 255' # Keep unique mappers
        + ' --threads ' + str(args_dict['threads'])
        + ' -o ' + str(args_dict['alignments']) + str(output) + '_Aligned.unique.bam'
        + ' ' + str(args_dict['alignments']) + str(output) + '_Aligned.out.bam'
        + str(args_dict['log']))

    # Sort SAM file
    os.system(
        'samtools sort'
        + ' --threads ' + str(args_dict['threads'])
        + ' -o ' + str(args_dict['alignments']) + str(output) + '_Aligned.sort.bam'
        + ' ' + str(args_dict['alignments']) + str(output) + '_Aligned.unique.bam'
        + str(args_dict['log']))

    # Index BAM file
    os.system(
        'samtools index'
        + ' ' + str(args_dict['alignments']) + str(output) + '_Aligned.sort.bam'
        + str(args_dict['log']))

    # Use sorted BAM file to find any duplicate reads
    os.system(
        'samtools markdup'
        + ' ' + str(args_dict['alignments']) + str(output) + '_Aligned.sort.bam' # Input BAM
        + ' ' + str(args_dict['alignments']) + str(output) + '_dedupMarked.bam' # Output BAM
        + ' -s' # Print some basic stats
        + str(args_dict['log']))
    os.system(
        'samtools index'
        + ' ' + str(args_dict['alignments']) + str(output) + '_dedupMarked.bam'
        + str(args_dict['log']))

    # Create sorted BAM file with duplicates removed
    os.system(
        'samtools markdup'
        + ' ' + str(args_dict['alignments']) + str(output) + '_Aligned.sort.bam' # Input BAM
        + ' ' + str(args_dict['alignments']) + str(output) + '_dedupRemoved.bam' # Output BAM
        + ' -s' # Print some basic stats to STDOUT
        + ' -r' # Remove duplicate reads
        + str(args_dict['log']))
    os.system(
        'samtools index'
        + ' ' + str(args_dict['alignments']) + str(output) + '_dedupRemoved.bam'
        + str(args_dict['log']))

"""Remove all intermediate alignment files and references after alignment is complete"""
def remove_intermediates(
        args_dict):

    os.system(
        "find"
        + " " + str(args_dict['alignments'])
        + " ! -name '*_Aligned.sort.bam'"
        + " ! -name '*_Aligned.sort.bam.bai'"
        + " ! -name '*_dedupMarked.bam'"
        + " ! -name '*_dedupRemoved.bam'"
        + " ! -name '*_dedupRemoved.bam.bai'"
        + " ! -name '*_Log.final.out'"
        + " -maxdepth 1 -type f -delete" # Only keep files matching pattern
        + str(args_dict['log']))

def clean_reference_directory(
        args_dict):

    os.system(
        'rm -r'
        + ' ' + str(args_dict['intermediate_references'])
        + str(args_dict['log']))

def align(
    args_dict,
    output,
    file):

    if 'mask' in args_dict and args_dict['mask'] == True:
        file = masking_star(
            file,
            output,
            args_dict)

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

    # Create BAM file with only unique hits, mark duplicates, index
    alignment_process(
        output,
        args_dict)

    # Clean up the output
    remove_intermediates(
        args_dict)

"""Single-end RNA-seq pipeline"""
def se_align(
        args):

    file, args_dict = args[0], args[1]

    output = str(file[8:-6]) # Get output file name before adding path to file name(s)
    file = str(args_dict['input']) + str(file)

    align(args_dict, output, file)

"""Paired-end RNA-seq pipeline"""
def pe_align(
        args):

    file1, file2, args_dict = args[0], args[1], args[2]

    # STAR first pass
    output = str(file1[8:-7]) # Get output file name before adding path to file name(s)
    file = str(args_dict['input']) + str(file1) + ' ' + str(args_dict['input']) + str(file2)

    align(args_dict, output, file)

"""Manage single-end RNA-seq pipeline"""
def run_seRNAseq(
        args_dict):

    try:
        # Add output directories
        args_dict = add_directory(
            args_dict,
            'output',
            'alignments')
        args_dict = add_directory(
            args_dict,
            'alignments',
            'intermediate_references')

        args_dict['fasta_list'] = get_fasta(args_dict['reference'])

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

    except:
        raise Exception('Single-end alignment failed')

"""Manage paired-end RNA-seq pipeline"""
def run_peRNAseq(
        args_dict):

    try:
        # Add output directories
        args_dict = add_directory(
            args_dict,
            'output',
            'alignments')
        args_dict = add_directory(
            args_dict,
            'alignments',
            'intermediate_references')

        args_dict['fasta_list'] = get_fasta(args_dict['reference'])

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

    except:
        raise Exception('Paired-end alignment failed')
