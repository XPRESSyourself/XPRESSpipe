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
import math
import concurrent.futures
from multiprocessing import cpu_count

"""Determine number of processors to use"""
def get_cores(
    args_dict, mod_workers):

    if 'max_processors' in args_dict and args_dict['max_processors'] != None:
        cores = args_dict['max_processors']
    else:
        cores = cpu_count() #Number of CPU cores on your system

    if mod_workers == True:
        workers = cores
    else:
        workers = 1

    return cores, workers

"""Run function and files on pools"""
def run_pools(func, args_iter, args_dict):

    with concurrent.futures.ProcessPoolExecutor(max_workers=args_dict['workers']) as executor:
        for file in zip(args_iter, executor.map(func, args_iter)):
            print(file, "has been processed.")

"""Parallelize function on list of files"""
def parallelize(func, file_list, args_dict, mod_workers=False):

    args_iter = ([file, args_dict] for file in file_list)

    args_dict['threads'], args_dict['workers'] = get_cores(args_dict, mod_workers)

    run_pools(func, args_iter, args_dict)

"""Parallelize function on list of files for PE data"""
def parallelize_pe(func, file_list, args_dict, mod_workers=False):

    # Pair files for paired-end processing
    c1 = 0
    args_iter = []
    for c in range(int(len(file_list)/2)):
        c2 = c1 + 1
        args_iter.append([file_list[c1], file_list[c2], args_dict])
        c1 += 2

    args_iter = ([x[0], x[1], x[2]] for x in args_iter)

    args_dict['threads'], args_dict['workers'] = get_cores(args_dict, mod_workers)

    run_pools(func, args_iter, args_dict)
