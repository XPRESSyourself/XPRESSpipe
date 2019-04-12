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

"""Read in indexed BAM file to Pandas dataframe"""
def read_BAM(
    bam_file, # Must be indexed before reading in
    threads=1):

    # Read in BAM file
    os.system('samtools view'
          + ' --threads ' + str(threads)
          + ' ' + str(bam_file)
          + ' -o ' + str(bam_file)[:-4] + '.tmp')

    #read sorted, unique only sam files, get most abundant length and make new file
    sam = pd.read_csv(
            str(bam_file)[:-4] + '.tmp',
            sep = '\t',
            header = None,
            usecols = list(range(0, 16)),
            low_memory = False)

    sam[2] = sam[2].astype(str)

    os.system('rm '
              + str(bam_file)[:-4] + '.tmp')

    return sam
