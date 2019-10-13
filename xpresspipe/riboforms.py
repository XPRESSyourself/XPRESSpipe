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

"""IMPORT DEPENDENCIES
"""
import os
import sys
import pandas as pd

"""Description:
Purpose of this module is to perform a Ribomap-like isoform quantification for ribosome profiling reads

Prereqs:
Run both HTSeq and Cufflinks (make this an option for quant by providing \"both\")
Quantification should have been performed using a GTF that was truncated and protein-coding-only

Inputs:
counts_table
abundance_table
meta_table -> should include two columns:
    index: sample names as in tables
    1. type

Outputs:
normalized_table

Protocol:
- Read in samples prefixes to list
- Parse through list and grab ribo count column and rna isoform abundance table for each
- For each gene in isoforms, get relative abundance of each transcript as dictionary
    d = {
        transcript_1: 0.9,
        transcript_2: 0.05,
        transcript_2: 0.05
    }
- Append transcript isoform table to master dataframe
- Calculate for each transcript gene counts * relative abundance and add to master dataframe
- Output table  
"""
