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

chromosome_position = 2
leftCoordinate_position = 3
read_position = 9

"""EXTERNAL DEPENDENCIES"""
# samtools

"""Read in indexed BAM file to Pandas dataframe"""
def read_bam(
    file,
    keep_unmapped = False,
    threads = 1):

    # Read in BAM file
    if str(file).lower()[-4:] == '.bam':
        os.system(
            'samtools view'
            + ' --threads ' + str(threads)
            + ' ' + str(file)
            + ' -o ' + str(file).rsplit('.',1)[0] + '.sam')

        #read sorted, unique only sam files, get most abundant length and make new file
        bam = pd.read_csv(
            str(file).rsplit('.',1)[0] + '.sam',
            sep = '\t',
            header = None,
            usecols = [chromosome_position, leftCoordinate_position, read_position],
            low_memory = False)

        os.system(
            'rm'
            ' ' + str(file).rsplit('.',1)[0] + '.sam')

    elif str(file).lower()[-4:] == '.sam':
        bam = pd.read_csv(
            str(file).rsplit('.',1)[0] + '.sam',
            sep = '\t',
            header = None,
            usecols = [chromosome_position, leftCoordinate_position, read_position],
            low_memory = False)

    else:
        raise Exception('A SAM/BAM format file was not provided.')

    bam[chromosome_position] = bam[chromosome_position].astype(str) # Make chromosome info consisitent with reference file
    if keep_unmapped == False:
        bam = bam.loc[bam[chromosome_position] != '*']

    return bam

"""Randomly sample n reads from BAM file"""
def bam_sample(
        bam,
        number):

    if bam.shape[0] < int(number):
        raise Exception('Not enough alignments to sample')

    return bam.sample(n = int(number))
