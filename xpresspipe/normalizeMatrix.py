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
from xpresstools import batch_normalize, rpm, tpm, r_fpkm, log_scale

"""Run normalization of count dataframe"""
def run_normalization(
        args_dict,
        sep='\t'):

    # Run sample normalization
    if 'method' in args_dict and args_dict['method'] != None:

        df = pd.read_csv(
            str(args_dict['input']),
            sep = sep,
            header = 0,
            index_col = 0,
            comment = '#',
            low_memory = False)

        # RPM normalization
        if args_dict['method'].upper() == 'RPM':
            type = 'rpm'
            df = rpm(df)
            df.to_csv(
                str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv',
                sep = '\t')

        # RPKM or FPKM normalization
        elif args_dict['method'].upper() == 'RPKM' \
        or args_dict['method'].upper() == 'FPKM':
            if args_dict['gtf'] == None:
                raise Exception('A GTF reference file is required for RPKM and FPKM normalization')
            type = 'r_fpkm'
            df = r_fpkm(
                df,
                args_dict['gtf'])
            df.to_csv(
                str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv',
                sep = '\t')
        elif args_dict['method'].upper() == 'TPM':
            if args_dict['gtf'] == None:
                raise Exception('A GTF reference file is required for RPKM and FPKM normalization')
            type = 'tpm'
            df = tpm(
                df,
                args_dict['gtf'])
            df.to_csv(
                str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv',
                sep = '\t')

        # Log normalization
        elif args_dict['method'].upper() == 'LOG':
            type = 'log'
            df = log_scale(
                df,
                log_base = 10)
            df.to_csv(
                str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.tsv',
                sep = '\t')
        else:
            raise Exception('Unknown \"method\" argument provided')

    # Run in batch normalization
        if 'batch' in args_dict and args_dict['batch'] != None:
            batch_normalize(
                str(args_dict['input'][:-4]) + '_' + str(type) + 'Normalized.csv',
                str(args_dict['batch']))
    else:
        if 'batch' in args_dict and args_dict['batch'] != None:
            batch_normalize(
                str(args_dict['input']),
                str(args_dict['batch']))
