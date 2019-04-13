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

"""IMPORT DEPENDENCIES"""
import os
import sys
import pandas as pd
import numpy as np
import math

"""EXTERNAL DEPENDENCIES"""
# samtools

"""Read in indexed BAM file to Pandas dataframe"""
def read_bam(file, threads = 1):

    # Read in BAM file
    os.system('samtools view'
          + ' --threads ' + str(threads)
          + ' ' + str(file)
          + ' -o ' + str(file)[:-4] + '.tmp')

    #read sorted, unique only sam files, get most abundant length and make new file
    bam = pd.read_csv(
            str(file)[:-4] + '.tmp',
            sep = '\t',
            header = None,
            usecols = list(range(0, 16)),
            low_memory = False)

    bam[2] = bam[2].astype(str)

    os.system('rm'
              ' ' + str(file)[:-4] + '.tmp')

    return bam

""""""
def bam_sample(bam, number):

    if bam.shape[0] < int(number):
        raise Exception('Not enough alignments to sample')

    return bam.sample(n = int(number))

""""""
def mid_coordinates(bam):

    bam[16] = bam[3] + (bam[9].str.len() / 2).apply(np.floor).astype('int64')
    mid_coordinates = bam[[2,16]]
    return mid_coordinates.values.tolist()

"""P site coordinate for 28-30mers"""
def phased_coordinates(bam):

    return phased_coordinates.values.tolist()
