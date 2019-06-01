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
import csv
import warnings
import pandas as pd
pd.options.mode.chained_assignment = None
import multiprocessing # For debugging
from multiprocessing import cpu_count, Pool
import gc
from functools import partial
from Bio import SeqIO

"""IMPORT INTERNAL DEPENDENCIES"""
from .gtfTruncate import truncate_gtf
from .utils import get_files

"""
Parse  GTF dataframe for longest transcript per gene record and keep only those transcript records
Input: GTF formatted as pandas dataframe
"""

def check_stops(
    data,
    stops):

    for index, row in data.iterrows():

        for s in stops:
            if s > row[3] or s < row[4]:
                return 1
            else:
                pass

    return 0


def longest_transcripts(
        gtf):

    long_transcripts = []

    for index, row in gtf.iterrows():

        # Find next gene record
        if row[2] == 'gene':
            # Forward scan to next gene record
            n = 0
            gene_id_original = gtf.at[index + n, 8][(gtf.at[index + n, 8].find('gene_id \"') + 9):].split('\";')[0]
            gene_id_next = gene_id_original

            while gene_id_next == gene_id_original:

                if index + n + 1 > gtf.index[-1]:
                    break
                else:
                    n += 1
                    gene_id_next = gtf.at[index + n, 8][(gtf.at[index + n, 8].find('gene_id \"') + 9):].split('\";')[0]

            # Parse out current gene record
            gene_start_index = index
            if gtf.at[index + n, 2] == 'gene':
                gene_stop_index = index + n - 1
            else:
                gene_stop_index = index + n

            if gene_start_index < gene_stop_index:
                # Get gene record
                gtf_record = gtf.loc[gene_start_index:gene_stop_index]

                # Get coordinates for each transcript record
                transcript_list = []
                for index, row in gtf_record.iterrows():
                    if row[2] == 'transcript':
                        transcript_list.append(index)
                transcript_list.append(gene_stop_index)

                # Get Priority, Length, Transcript ID for each transcript record
                transcript_info = []
                for x in range(0, len(transcript_list)-1):
                    if transcript_list[x + 1] == gene_stop_index:
                        transcript_record = gtf.loc[transcript_list[x]:transcript_list[x + 1]]
                    else:
                        transcript_record = gtf.loc[transcript_list[x]:transcript_list[x + 1] - 1]

                    # Get transcript ID
                    transcript_id = transcript_record.iloc[0,8][(transcript_record.iloc[0, 8].find('transcript_id \"') + 15):].split('\";')[0]

                    # For each transcript record, get stop codon coordinates, CDS and
                    stops = []
                    for index, row in transcript_record.iterrows():
                        if row[2] == 'stop_codon':
                            stops.append(row[3])
                            stops.append(row[3] + 1)
                            stops.append(row[4])

                    # Get values of interest for each transcript record and check if stop in range
                    cds = transcript_record.loc[transcript_record[2] == 'CDS']
                    if check_stops(cds, stops) == 1:
                        cds = pd.DataFrame()

                    exon_ens_hav = transcript_record.loc[(transcript_record[2] == 'exon') & (transcript_record[1].isin(['ensembl', 'havana', 'ensembl_havana']))]
                    if check_stops(exon_ens_hav, stops) == 1:
                        exon_ens_hav = pd.DataFrame()

                    exon = transcript_record.loc[transcript_record[2] == 'exon']
                    if check_stops(exon, stops) == 1:
                        priority = 4
                    else:
                        priority = 3

                    # Set priority and length for each type
                    length = 0
                    if cds.empty == False:
                        priority = 1
                        for index, row in cds.iterrows():
                            length = length + (row[4] - row[3])
                    elif exon_ens_hav.empty == False:
                        priority = 2
                        for index, row in exon_ens_hav.iterrows():
                            length = length + (row[4] - row[3])
                    elif exon.empty == False:
                        for index, row in exon_ens_hav.iterrows():
                            length = length + (row[4] - row[3])
                    else:
                        raise Exception('Something went wrong')

                    transcript_info.append([priority, length, transcript_id])

                # Compare the different transcripts for highest priority and longest transcript
                priorities = []
                for x in transcript_info:
                    priorities.append(x[0])
                priority_max = max(priorities)

                lengths = []
                for x in transcript_info:
                    if x[0] == priority_max:
                        lengths.append(x)
                lengths_max = max(lengths)
                print(lengths_max)
                for x in transcript_info:
                    if x[1] == lengths_max:
                        print(x)
                        long_transcripts.append(gtf.loc[gtf[8].str.contains(x[1])])
    print(long_transcripts)

    #gtf = None # Garbage management
    gc.collect()

    if len(long_transcripts) > 0:
        gtf_longest = pd.concat(long_transcripts)
        gtf_longest = gtf_longest.reset_index(drop=True)

        long_transcripts = None
        gc.collect()

        return gtf_longest

    else:
        return

"""Only keep GTF records that are annotated as protein_coding"""
def protein_gtf(
        gtf):

    # Take only records that are annotated as 'protein coding'
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]
    gtf_coding = gtf_coding.reset_index(drop=True)

    gtf = None # Garbage management
    gc.collect()

    return gtf_coding

"""
Parse GTF down to chunks per cores for multiprocessing
Requires unmodified GTF with gene records intact
"""
def get_chunks(
        gtf,
        threads=None):

    # Determine number of chunks to create based on indicated or available cpus
    cores = cpu_count() # Number of CPU cores on your system
    if threads == None or threads >= cores:
        pass
    elif threads < 1:
        warnings.warn('Indicated less than 1 CPU, setting to 1')
        cores = 1
    else:
        cores = int(threads)

    # Get number of times 'gene' is references in column 2 of GTF and if cores > #genes, limit
    gene_instances = gtf[gtf[2] == 'gene'][2].shape[0]
    if gene_instances == 0:
        raise Exception('No gene records found in GTF')

    if cores > gene_instances:
        cores = gene_instances

    start = 0 # Get first start coordinate for chunk
    batch = round(len(gtf.index) / cores) # Approx. number of samples in a chunk

    # Get chunking indices
    chunks = [] # Initialize chunk storage
    for y in range(cores):

        # If the last chunk, get the remainder of the GTF dataframe
        if y == cores - 1:
            new_chunk = gtf.loc[start:]

        else:
            end = start + batch  # Set tentative end of next chunk

            if end > len(gtf.index) - 1: # If end of dataframe, end there
                end = len(gtf.index) - 1

            else: # Go to tentative end of next chunk and search until previous gene record
                gtf_remainder = gtf.iloc[start:end]

                n = -1 # Start at current record
                gene_id_original = gtf_remainder.at[end + n, 8][(gtf_remainder.at[end + n, 8].find('gene_id \"') + 9):].split('\";')[0]
                gene_id_next = gene_id_original

                while gene_id_next == gene_id_original:
                    n -= 1 # Take another step back until a new gene_id is found
                    gene_id_next = gtf_remainder.at[end + n, 8][(gtf_remainder.at[end + n, 8].find('gene_id \"') + 9):].split('\";')[0]

            # Parse out current chunk
            new_chunk = gtf.loc[start:end + n]
            start = end + n + 1 # End coordinate for last chunk to start with next

        if new_chunk.empty == False:
            chunks.append(new_chunk)

    print('Dataframe split into ' + str(len(chunks)) + ' chunks for parallelization')
    gtf = None # Garbage management
    gc.collect()

    return chunks

"""Run a given function on chunks"""
def run_chunks(
        func,
        chunks,
        target_message=''):

    if target_message == '':
        print('Parsing record')
    else:
        print('Parsing record for ' + str(target_message))

    chunks = [x for x in chunks if x is not None] # Remove any empty dataframes
    cores = len(chunks) # Modify worker numbers
    pool = Pool(cores) # Initialize workers
    chunks = pool.map(func, chunks) # Run function on chunks
    pool.close()
    pool.join()
    gc.collect()

    return chunks

"""Run all GTF-editing functions"""
def edit_gtf(
    gtf, # Dataframe of file path and name to GTF reference
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=True,
    _5prime=45, # If no 5' truncation desired, set to 0
    _3prime=15, # If no 3' truncation desired, set to 0
    output=True, # True will output all intermediates, not possible if inputting a GTF as pandas dataframe
    threads=None): # Give int for core threshold if desired

    # Import GTF reference file
    if isinstance(gtf, pd.DataFrame) and len(gtf.columns) == 9:
        output = False # Turn off intermediates output
        file_name = None
    elif str(gtf).endswith('.gtf'):
        file_name = str(gtf[:-4]) + '_' # Get rid of GTF extension for now
        gtf = pd.read_csv(
            str(gtf),
            sep = '\t',
            header = None,
            comment = '#',
            low_memory = False)
    else:
        raise Exception('Error: A GTF-formatted file or dataframe was not provided')

    # Get chunks
    chunks = get_chunks(
        gtf,
        threads = threads)

    gtf = None # Garbage management
    gc.collect()

    # Run GTF modifications
    # Parse each gene record for longest transcript
    if longest_transcript == True:
        chunks = run_chunks(
            longest_transcripts,
            chunks,
            target_message = 'longest transcripts')

        if output == True:
            file_name = str(file_name) + 'L'

    # Get only protein coding annotated records
    if protein_coding == True:
        chunks = run_chunks(
            protein_gtf,
            chunks,
            target_message = 'protein coding genes')

        if output == True:
            file_name = str(file_name) + 'C'

    # Truncate by unique transcript
    # If file has not been parsed for longest transcript per gene, will truncate each isoform
    if truncate_reference == True:
        func = partial(
            truncate_gtf,
            _5prime = _5prime,
            _3prime = _3prime)
        chunks = run_chunks(
            func,
            chunks,
            target_message = 'truncation')

        if output == True:
            file_name = str(file_name) + 'T'

    # Merge final GTF from chunks and output
    if len(chunks) > 0:
        gtf = pd.concat(chunks)
        gtf = gtf.reset_index(drop=True)

        chunks = None # Garbage management
        gc.collect()

        if output == True:
            gtf.to_csv(
                str(file_name) + '.gtf',
                sep = '\t',
                header = None,
                index = False,
                quoting = csv.QUOTE_NONE)
        else:
            return gtf

    else:
        raise Warning('0 chunks of the original file remain')
