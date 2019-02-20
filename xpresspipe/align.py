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
from .utils import get_files, unzip_files, add_directory
from .parallel import parallelize, parallelize_pe, compute_cores_files

"""
DESCRIPTION: Remove all intermediate alignment files and references after alignment is complete
"""
def remove_intermediates(args_dict):

    os.system('rm -r ' + str(args_dict['intermediate_references']))
    os.system('rm !(' + str(args_dict['input']) + '*_final.sam|' + str(args_dict['input']) + '*_final_Log.final.out)') #removes all but final alignment file

"""
DESCRIPTION:
"""
def se_align(args):

    file, args_dict = args[0], args[1]

    #STAR first pass
    os.system('STAR --genomeDir ' + str(args_dict['reference']) + ' --readFilesIn ' + str(args_dict['input']) + str(file) + ' --runThreadN ' + str(args_dict['threads']) + ' --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None --outFileNamePrefix ' + str(args_dict['alignments']) + str(file[:-6]) + '_')

    #STAR intermediate reference building
    os.system('mkdir ' + str(args_dict['intermediate_references']) + str(file[:-6]))
    os.system('STAR --runMode genomeGenerate --genomeDir ' + str(args_dict['intermediate_references']) + str(file[:-6]) + ' --genomeFastaFiles ' + str(args_dict['intermediate_references']) + 'genome.fa --sjdbOverhang 100 --runThreadN ' + str(args_dict['threads']) + ' --sjdbFileChrStartEnd ' + str(args_dict['alignments']) + str(file[:-6]) = '_SJ.out.tab')

    #STAR second pass
    os.system('STAR --genomeDir ' + str(args_dict['intermediate_references']) + str(file[:-6]) + ' --readFilesIn ' + str(args_dict['input']) + str(file) + ' --runThreadN ' + str(args_dict['threads']) + ' --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype SAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outFileNamePrefix ' + str(args_dict['alignments']) + str(file[:-6]) + '_final_')

    #Create sam file with only unique hits
    os.system('samtools sort ' + str(args_dict['alignments']) + str(file[:-6]) + '_final_Aligned.sam -o ' + str(args_dict['alignments']) + str(file[:-6]) + '_sorted.sam')
    os.system('samtools view -q 255 ' + str(args_dict['alignments']) + str(file[:-6]) + '_sorted.sam > ' + str(args_dict['alignments']) + str(file[:-6]) + '_final.sam')

"""
DESCRIPTION:
"""
def pe_align(args):

    file1, file2, args_dict = args[0], args[1], args[2]

    #STAR first pass
    os.system('STAR --genomeDir ' + str(args_dict['reference']) + ' --readFilesIn ' + str(args_dict['input']) + str(file1) + ' ' + str(args_dict['input']) + str(file2) + ' --runThreadN ' + str(args_dict['threads']) + ' --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None --outFileNamePrefix ' + str(args_dict['alignments']) + str(file1[:-6]) + '_')

    #STAR intermediate reference building
    os.system('mkdir ' + str(args_dict['intermediate_references']) + str(file1[:-6]))
    os.system('STAR --runMode genomeGenerate --genomeDir ' + str(args_dict['intermediate_references']) + str(file1[:-6]) + ' --genomeFastaFiles ' + str(args_dict['intermediate_references']) + 'genome.fa --sjdbOverhang 100 --runThreadN ' + str(args_dict['threads']) + ' --sjdbFileChrStartEnd ' + str(args_dict['alignments']) + str(file1[:-6]) = '_SJ.out.tab')

    #STAR second pass
    os.system('STAR --genomeDir ' + str(args_dict['intermediate_references']) + str(file1[:-6]) + ' --readFilesIn ' + str(args_dict['input']) + str(file1) + ' ' + str(args_dict['input']) + str(file2) + ' --runThreadN ' + str(args_dict['threads']) + ' --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype SAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outFileNamePrefix ' + str(args_dict['alignments']) + str(file1[:-6]) + '_final_')

    #Create sam file with only unique hits
    os.system('samtools sort ' + str(args_dict['alignments']) + str(file1[:-6]) + '_final_Aligned.sam -o ' + str(args_dict['alignments']) + str(file1[:-6]) + '_sorted.sam')
    os.system('samtools view -q 255 ' + str(args_dict['alignments']) + str(file1[:-6]) + '_sorted.sam > ' + str(args_dict['alignments']) + str(file1[:-6]) + '_final.sam')

    #Do paired reads need to be collated? -- downstream stuff assumes yes


"""
DESCRIPTION: Align single-end Illumina RNAseq reads
"""
def run_seRNAseq(args_dict):

    try:
        #Add output directories
        args_dict = add_directory(args_dict, 'output', 'alignments')
        args_dict = add_directory(args_dict, 'alignments', 'intermediate_references')

        #Unzip files
        unzip_files(args_dict['input'])

        #Get list of files to align based on acceptable file types
        files = get_files(args_dict['input'], ['.fastq','.fq','.txt'])

        #Align single-end RNAseq reads
        parallize(se_align, files, args_dict)
        remove_intermediates(args_dict)

        return args_dict

    except:
        raise Exception('Paired-end alignment failed')

"""
DESCRIPTION: Align paired-end Illumina RNAseq reads
"""
def run_peRNAseq(args_dict):

    try:
        #Add output directories
        args_dict = add_directory(args_dict, 'output', 'alignments')
        args_dict = add_directory(args_dict, 'alignments', 'intermediate_references')

        #Unzip files
        unzip_files(args_dict['input'])

        #Get list of files to align based on acceptable file types
        files = get_files(args_dict['input'], ['.fastq','.fq','.txt'])

        if len(files) % 2 != 0:
            raise Exception('An uneven number of paired-end files were specified in the input directory')
        else:
            #Align paired-end RNAseq reads
            parallelize_pe(pe_align, files, args_dict)
            remove_intermediates(args_dict)

        return args_dict

    except:
        raise Exception('Paired-end alignment failed')
