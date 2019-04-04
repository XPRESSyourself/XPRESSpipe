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

"""IMPORT DEPENDENCIES"""
import os
import sys

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import get_files, unzip_files, add_directory, get_fasta
from .parallel import parallelize, parallelize_pe

"""Run first pass STAR alignment to map splice junctions"""
def first_pass_star(
    file, output, args_dict):

    os.system('STAR' \
        + ' --runThreadN ' + str(args_dict['threads']) \ # Argument to specify number of threads to use for processing
        + ' --genomeDir ' + str(args_dict['reference']) \ # Argument for specifying STAR reference directory
        + 'genome --readFilesIn ' + str(file) \ # Argument to dictate directory where pre-processed read files are located
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_' \ # Argument to dictate output directory
        + ' --seedSearchStartLmax ' + str(args_dict['seedSearchStartLmax']) \ #
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) \ # Read overhand amount to allow for splice mapping (should be same used in curation of reference)
        + ' --outFilterMultimapScoreRange 1' \
        + ' --outFilterMultimapNmax 20' \
        + ' --outFilterMismatchNmax 10' \
        + ' --alignIntronMax 500000' \
        + ' --alignMatesGapMax 1000000' \
        + ' --sjdbScore 2' \
        + ' --alignSJDBoverhangMin 1' \
        + ' --genomeLoad NoSharedMemory' \
        + ' --readFilesCommand cat' \
        + ' --outFilterMatchNminOverLread 0.33' \
        + ' --outFilterScoreMinOverLread 0.33' \
        + ' --outSAMstrandField intronMotif' \
        + ' --outSAMtype None' \
        + ' --outSAMmode None' \
        + str(args_dict['log']) \ # Record log output (must go last in command)
        )

"""Build intermediate STAR alignment reference using splice junction annotations from first pass"""
def build_star_splice_junction_intermediate(
    output, args_dict):

    os.system('mkdir' \
        + ' ' \
        + str(args_dict['intermediate_references']) + str(output) \ # Make directory for currently processed file
        + str(args_dict['log']) \ # Record log output (must go last in command)
        )
    os.system('STAR' \
        + ' --runMode genomeGenerate' \
        + ' --runThreadN ' + str(args_dict['threads']) \
        + ' --sjdbFileChrStartEnd ' + str(args_dict['alignments']) + str(output) + '_SJ.out.tab' \ # Splice junction annotation file output in first pass alignment
        + ' --genomeFastaFiles ' + str(args_dict['fasta_list']) \ # Input chromosomal fasta files for reference building
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output) \ # Location for output revised reference
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) \ # Read overhand amount to allow for splice mapping (should be same used in curation of reference)
        + str(args_dict['log']) \ # Record log output (must go last in command)
        )

"""Run second pass STAR alignment to map reads splice-aware"""
def second_pass_star(
    file, output, args_dict):

    os.system('STAR' \
        + ' --runThreadN ' + str(args_dict['threads']) \
        + ' --genomeDir ' + str(args_dict['intermediate_references']) + str(output) \
        + ' --readFilesIn ' + str(file) \
        + ' --outFileNamePrefix ' + str(args_dict['alignments']) + str(output) + '_final_' \
        + ' --seedSearchStartLmax ' + str(args_dict['seedSearchStartLmax']) \
        + ' --outFilterMultimapScoreRange 1' \
        + ' --outFilterMultimapNmax 20' \
        + ' --outFilterMismatchNmax 10' \
        + ' --alignIntronMax 500000' \
        + ' --alignMatesGapMax 1000000' \
        + ' --sjdbScore 2' \
        + ' --alignSJDBoverhangMin 1' \
        + ' --genomeLoad NoSharedMemory' \
        + ' --limitBAMsortRAM 0' \
        + ' --readFilesCommand cat' \
        + ' --outFilterMatchNminOverLread 0.33' \
        + ' --outFilterScoreMinOverLread 0.33' \
        + ' --sjdbOverhang ' + str(args_dict['sjdbOverhang']) \
        + ' --outSAMstrandField intronMotif' \
        + ' --outSAMattributes NH HI NM MD AS XS' \
        + ' --outSAMunmapped Within' \
        + ' --outSAMtype SAM' \
        + ' --outSAMheaderHD @HD VN:1.4' \
        + str(args_dict['log']) \
        )

"""Sort reads per file by chromosome position and keep only unique mappers"""
def alignment_sort(
    output, args_dict):

    # Sort by coordinate
    os.system('samtools sort' \
        + ' --threads ' + str(args_dict['threads']) \
        + ' ' + str(args_dict['alignments']) + str(output) + '_final_Aligned.out.sam' \
        + ' -o ' + str(args_dict['alignments']) + str(output) + '_sorted.sam' \
        + str(args_dict['log']) \
        )

    # only take unique mappers (q = 255)
    os.system('samtools view' \
        + ' -h -q 255' \
        + ' --threads ' + str(args_dict['threads']) \
        + ' ' + str(args_dict['alignments']) + str(output) + '_sorted.sam' \
        + ' > ' + str(args_dict['alignments']) + str(output) + '_final.sam' \
        )

    # make bam file
    os.system('samtools view' \
        + ' -S -b' \
        + ' --threads ' + str(args_dict['threads']) \
        + ' ' + str(args_dict['alignments']) + str(output) + '_final.sam' \
        + ' > ' + str(args_dict['alignments']) + str(output) + '_final.bam'
        )

"""Remove all intermediate alignment files and references after alignment is complete"""
def remove_intermediates(args_dict):

    os.system('find' \
        + ' ' + str(args_dict['alignments']) + ' ! -name *_final.bam ! -name *_final.sam ! -name *_final_Log.final.out -maxdepth 1 -type f -delete' \ # Only keep files matching pattern
        + str(args_dict['log'])
        )

def clean_directories(args_dict):

    os.system('rm -r' \
        + ' ' + str(args_dict['intermediate_references']) \
        + str(args_dict['log'])
        )

"""Single-end RNA-seq pipeline"""
def se_align(args):

    file, args_dict = args[0], args[1]

    #STAR first pass
    output = str(file[:-6]) #get output file name before adding path to file name(s)
    file = str(args_dict['input']) + str(file)
    first_pass_star(file, output, args_dict)

    #STAR intermediate reference building
    build_star_splice_junction_intermediate(output, args_dict)

    #STAR second pass
    second_pass_star(file, output, args_dict)

    #Create sam file with only unique hits
    alignment_sort(output, args_dict)

    remove_intermediates(args_dict)

"""Paired-end RNA-seq pipeline"""
def pe_align(args):

    file1, file2, args_dict = args[0], args[1], args[2]

    #STAR first pass
    output = str(file1[:-7]) #get output file name before adding path to file name(s)
    file = str(args_dict['input']) + str(file1) + ' ' + str(args_dict['input']) + str(file2)

    first_pass_star(file, output, args_dict)

    #STAR intermediate reference building
    build_star_splice_junction_intermediate(output, args_dict)

    #STAR second pass
    second_pass_star(file, output, args_dict)

    #Create sam file with only unique hits
    alignment_sort(output, args_dict)

    remove_intermediates(args_dict)

"""Manage single-end RNA-seq pipeline"""
def run_seRNAseq(args_dict):

    try:
        #Add output directories
        args_dict = add_directory(args_dict, 'output', 'alignments')
        args_dict = add_directory(args_dict, 'alignments', 'intermediate_references')

        args_dict['fasta_list'] = get_fasta(args_dict['reference'])

        #Unzip files
        unzip_files(args_dict['input'])

        #Get list of files to align based on acceptable file types
        files = get_files(args_dict['input'], ['.fastq','.fq','.txt'])

        #Align single-end RNAseq reads
        parallelize(se_align, files, args_dict)
        clean_directories(args_dict)

        return args_dict

    except:
        raise Exception('Single-end alignment failed')

"""Manage paired-end RNA-seq pipeline"""
def run_peRNAseq(args_dict):

    try:
        #Add output directories
        args_dict = add_directory(args_dict, 'output', 'alignments')
        args_dict = add_directory(args_dict, 'alignments', 'intermediate_references')

        args_dict['fasta_list'] = get_fasta(args_dict['reference'])

        #Unzip files
        unzip_files(args_dict['input'])

        #Get list of files to align based on acceptable file types
        files = get_files(args_dict['input'], ['.fastq','.fq','.txt'])

        if len(files) % 2 != 0:
            raise Exception('An uneven number of paired-end files were specified in the input directory')
        else:
            #Align paired-end RNAseq reads
            parallelize_pe(pe_align, files, args_dict)
            clean_directories(args_dict)

        return args_dict

    except:
        raise Exception('Paired-end alignment failed')
