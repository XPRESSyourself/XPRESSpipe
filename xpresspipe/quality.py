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

"""IMPORT INTERNAL DEPENDENCIES"""
from .utils import add_directory, get_files
from .parallel import parallelize

"""Generate read distribution profiles"""
def run_fastqc(
        args):

    file, args_dict = args[0], args[1]

    os.system(
        'fastqc'
        + ' -q ' + str(args_dict['input']) + str(file)
        + ' -o ' + str(args_dict['fastqc'])
        + str(args_dict['log']))

"""Parellel run of FastQC"""
def get_fastqc(args_dict):

    args_dict = add_directory(
        args_dict,
        'output',
        'fastqc')

    files = get_files(
        args_dict['input'],
        ['.fastq', '.fq', '.txt'])

    # Perform fastqc on each file and unzip output
    parallelize(
        run_fastqc,
        files,
        args_dict,
        mod_workers = True)

"""Create MultiQC processing summary from all files in args_dict output"""
def get_multiqc_summary(
        args_dict):

    os.system(
        'multiqc'
        + ' ' + str(args_dict['output'])
        + ' -i ' + str(args_dict['experiment'])
        + ' -o ' + args_dict['output']
        + str(args_dict['log']))
