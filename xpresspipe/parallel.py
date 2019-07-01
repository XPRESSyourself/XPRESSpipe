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
import math
import concurrent.futures
from multiprocessing import cpu_count
import psutil
import resource

"""
Func: Threshold number of workers if available RAM is insufficient with number of workers and file sizes
@param args_dict: Argument dictionary
@param file_list: List of file names used to check for max file size
- Get max file size in list and total memory available
- If it seems like the expansion of memory will kill the process, modify the number of workers allowed at a time to ensure RAM doesn't get overused at once
- Assume incoming BAM files will expand a lot more than a normal, non-binary file
"""
def threshold_ram(
    args_dict,
    file_list):

    total = psutil.virtual_memory()[1] # Get available memory

    file_sizes = [] # Get max file size
    for file in file_list:
        file_sizes.append(os.path.getsize(str(args_dict['input']) + str(file)))

    _max = max(file_sizes)

    if file[-4:] == '.bam': # Assume binary files will expand by factor of 4 for decompression and additional data storage used in process
        factor = 4
    else:
        factor = 1.5

    threshold_workers = math.floor(total / (_max * factor)) # Set threshold based on max file size in set
    if threshold_workers < 1:
        threshold_workers = 1

    if threshold_workers > cpu_count():
        threshold_workers = cpu_count()

    if threshold_workers < args_dict['workers']: # Modify if set # of workers is greater than memory threshold
        threshold_threads = math.floor(args_dict['threads'] / threshold_workers)
        print('Resetting parallelization specs based on max file size to be processed:\nMax number of workers: ' + str(threshold_workers) + '\nNumber of threads per worker (where available): ' + str(threshold_threads))
        return threshold_threads, threshold_workers
    else:
        return args_dict['threads'], args_dict['workers']

"""
Func: Determine number of processors to use
@param args_dict: Argument dictionary
@param mod_workers: Call to allow number of workers be equal to number of processors, else process one file at a time with all available processors
- Check number given as max processors and use that if not None
- If None specified (no user input), set cores equal to number available on system
- Determine number of workers to use per job based on user input
    - If modified, workers equal number of cores
    - If not modified, workers equal 1, so one worker is using all available cores
"""
def get_cores(
    args_dict,
    mod_workers):

    if 'max_processors' in args_dict and args_dict['max_processors'] != None:
        cores = args_dict['max_processors']
    else:
        cores = cpu_count() #Number of CPU cores on your system

    if mod_workers == True:
        workers = cores
    else:
        workers = 1

    return cores, workers

"""
Func: Run function and files on pools
@param func: function name to be executed on every object passed to the pool
@param args_iter: List of lists of file name and args_dict
@param args_dict: Argument dictionary
- Create batches of n args_iter objects. Each batch based on number of workers available at a given time
- Concurrently execute each file within a batch, clean process memory, pass in next batch, etc.
"""
def run_pools(
    func,
    args_iter,
    args_dict):

    pools = math.ceil(len(args_iter) / args_dict['workers'])

    it_list = []
    range_number = 0
    for x in range(pools):
        it_list.append([iter for iter in args_iter[range_number:range_number + args_dict['workers']]])
        range_number += args_dict['workers']

    batch_number = 1
    for batch in it_list:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args_dict['workers']) as executor:
            for file in zip(batch, executor.map(func, batch)):
                print(file[0][0], "has been processed.")
        print('Processing of batch {0} of {1} complete...'.format(batch_number, pools))
        batch_number += 1

"""
Func: Parallelize function on list of files
@param func: function name to be executed on every object passed to the pool
@param file_list: List of file names used to process
@param args_dict: Argument dictionary
@param mod_workers: Call to allow number of workers be equal to number of processors, else process one file at a time with all available processors
"""
def parallelize(
    func,
    file_list,
    args_dict,
    mod_workers=False):

    args_iter = [[file, args_dict] for file in file_list]

    # Get number of cores
    args_dict['threads'], args_dict['workers'] = get_cores(
        args_dict,
        mod_workers)

    # Check and apply RAM threshold if necessary
    if mod_workers == True:
        args_dict['threads'], args_dict['workers'] = threshold_ram(
            args_dict,
            file_list)

    run_pools(
        func,
        args_iter,
        args_dict)

"""
Func: Parallelize function on list of files for PE data
@param func: function name to be executed on every object passed to the pool
@param file_list: List of file names used to process
@param args_dict: Argument dictionary
@param mod_workers: Call to allow number of workers be equal to number of processors, else process one file at a time with all available processors
"""
def parallelize_pe(
    func,
    file_list,
    args_dict,
    mod_workers=False):

    # Pair files for paired-end processing
    c1 = 0
    args_iter = []
    for c in range(int(len(file_list)/2)):
        c2 = c1 + 1
        args_iter.append([file_list[c1], file_list[c2], args_dict])
        c1 += 2

    args_iter = [[x[0], x[1], x[2]] for x in args_iter]

    args_dict['threads'], args_dict['workers'] = get_cores(
        args_dict,
        mod_workers)

    # Check and apply RAM threshold if necessary
    if mod_workers == True:
        args_dict['threads'], args_dict['workers'] = threshold_ram(
            args_dict,
            file_list)

    run_pools(
        func,
        args_iter,
        args_dict)
