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

"""GLOBALS"""
gtf_chr_column = 1
gtf_type_column = 2
gtf_leftCoordinate_column = 3
gtf_rightCoordinate_column = 4
gtf_sign_column = 6
gtf_annotation_column = 8
parse_type = 'transcript_id \"'

"""
Parse  GTF dataframe for longest transcript per gene record and keep only those transcript records
Input: GTF formatted as pandas dataframe
"""
# Check if stop codon coordinates are within any CDS or exon space
def check_stops(
    data,
    stops):

    for index, row in data.iterrows():

        for s in stops:
            if s > row[gtf_leftCoordinate_column] and s < row[gtf_rightCoordinate_column]:
                return 1
            else:
                pass

    return 0

# Find longest trascript for each
# In case of tie breakers, first listed is taken
def longest_transcripts(
        gtf):

    long_transcripts = []

    for index, row in gtf.iterrows():

        # Find next gene record
        if row[gtf_type_column].lower() == 'gene':
            # Forward scan to next gene record
            n = 0
            gene_id_original = gtf.at[index + n, gtf_annotation_column][(gtf.at[index + n, gtf_annotation_column].find('gene_id \"') + 9):].split('\";')[0]
            gene_id_next = gene_id_original

            while gene_id_next == gene_id_original:

                if index + n + 1 > gtf.index[-1]:
                    break
                else:
                    n += 1
                    gene_id_next = gtf.at[index + n, gtf_annotation_column][(gtf.at[index + n, gtf_annotation_column].find('gene_id \"') + 9):].split('\";')[0]

            # Parse out current gene record
            gene_start_index = index
            if gtf.at[index + n, gtf_type_column].lower() == 'gene':
                gene_stop_index = index + n - 1
            else:
                gene_stop_index = index + n

            if gene_start_index < gene_stop_index:
                # Get gene record
                gtf_record = gtf.loc[gene_start_index:gene_stop_index]

                # Get coordinates for each transcript record
                transcript_list = []
                for index, row in gtf_record.iterrows():
                    if row[gtf_type_column].lower() == 'transcript':
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
                    transcript_id = transcript_record.iloc[0,gtf_annotation_column][(transcript_record.iloc[0, gtf_annotation_column].find('transcript_id \"') + 15):].split('\";')[0]

                    # For each transcript record, get stop codon coordinates, CDS and
                    stops = []
                    for index, row in transcript_record.iterrows():
                        if row[gtf_type_column].lower() == 'stop_codon':
                            stops.append(row[gtf_leftCoordinate_column])
                            stops.append(row[gtf_leftCoordinate_column] + 1)
                            stops.append(row[gtf_rightCoordinate_column])

                    # Get CCDS
                    #print(transcript_record)
                    ccds = transcript_record.loc[(transcript_record[gtf_type_column].str.upper() == 'CDS') & (transcript_record[gtf_annotation_column].str.contains('tag \"CCDS\"', case=False))]

                    if check_stops(ccds, stops) == 1:
                        ccds = pd.DataFrame()

                    trans_ens_hav = transcript_record.loc[(transcript_record[gtf_type_column].str.upper() == 'CDS') & (transcript_record[1].str.lower() == 'ensembl_havana')]
                    if check_stops(trans_ens_hav, stops) == 1:
                        trans_ens_hav = pd.DataFrame()

                    trans = transcript_record.loc[transcript_record[gtf_type_column].str.upper() == 'CDS']
                    if check_stops(trans, stops) == 1:
                        trans = pd.DataFrame()

                    exon_processed = transcript_record.loc[(transcript_record[gtf_type_column].str.lower() == 'exon') & (transcript_record[gtf_annotation_column].str.contains('transcript_biotype \"processed_transcript\"', case=False))]

                    exon = transcript_record.loc[transcript_record[gtf_type_column].str.lower() == 'exon']

                    # Set priority and length for each type
                    length = 0
                    if ccds.empty == False:
                        priority = 1
                        for index, row in ccds.iterrows():
                            length = length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1
                    elif trans_ens_hav.empty == False:
                        priority = 2
                        for index, row in trans_ens_hav.iterrows():
                            length = length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1
                    elif trans.empty == False:
                        priority = 3
                        for index, row in trans.iterrows():
                            length = length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1
                    elif exon_processed.empty == False:
                        priority = 4
                        for index, row in exon_processed.iterrows():
                            length = length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1
                    elif exon.empty == False:
                        priority = 5
                        for index, row in exon.iterrows():
                            length = length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1
                    else:
                        priority = 6
                        length = 0

                    # Get overall transcribed length for CCDS tie-breakers
                    transcript_length = 0
                    for index, row in exon.iterrows():
                        transcript_length = transcript_length + (row[gtf_rightCoordinate_column] - row[gtf_leftCoordinate_column]) + 1

                    # Make decision dataframe
                    transcript_info.append([priority, length, transcript_length, transcript_id])
                    transcript_data = pd.DataFrame.from_records(transcript_info)

                # Compare the different transcripts for highest priority and longest transcript (each successive step takes cohort meeting previous step)
                priority_max = min(transcript_data[0].tolist())
                coding_max = max(transcript_data.loc[transcript_data[0] == priority_max][1].tolist()) # From all records having the best priority
                transcript_max = max(transcript_data.loc[(transcript_data[0] == priority_max) & (transcript_data[1] == coding_max)][gtf_type_column].tolist()) # From all records with best priority and coding length

                # Keep first record that meets all above maxes in case of remaining tie-breaker
                transcript_keep = transcript_data.loc[(transcript_data[0] == priority_max) & (transcript_data[1] == coding_max) & (transcript_data[gtf_type_column] == transcript_max)][gtf_leftCoordinate_column].tolist()[0]
                long_transcripts.append(gtf.loc[gtf[gtf_annotation_column].str.contains(transcript_keep, case=False)])

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
    gtf_coding = gtf[gtf.iloc[:, gtf_annotation_column].str.contains('protein_coding', case=False) == True]
    gtf_coding = gtf_coding.reset_index(drop=True)

    gtf = None # Garbage management
    gc.collect()

    return gtf_coding

def get_ucsc_proteins(
        gtf,
        identifier='transcript_id \"',
        type='CDS'):

    gtf['gene'] = gtf[8].str.split(identifier).str[1]
    gtf['gene'] = gtf['gene'].str.split('\";').str[0]
    protein_list = set(gtf.loc[gtf[2] == type]['gene'].tolist())

    gtf_coding = gtf.copy()
    gtf_coding = gtf_coding[gtf_coding['gene'].isin(protein_list)]
    gtf_coding = gtf_coding.reset_index(drop=True)
    gtf_coding = gtf_coding.drop('gene', axis=1)

    return gtf_coding

"""
Parse GTF down to chunks per cores for multiprocessing
Requires unmodified GTF with gene records intact
"""
def get_chunks(
    gtf,
    threads=None,
    identifier='gene_id \"'):

    # Determine number of chunks to create based on indicated or available cpus
    cores = cpu_count() # Number of CPU cores on your system
    if threads == None or threads >= cores:
        pass
    elif threads < 1:
        print('Warning: Indicated less than 1 CPU, setting to 1')
        cores = 1
    else:
        cores = int(threads)

    # Get number of times 'gene' is references in column 2 of GTF and if cores > #genes, limit
    if len(gtf.columns.tolist()) == 9:
        gtf_check = gtf.copy()
        gtf_check['gene'] = gtf_check[8].str.split(identifier).str[1]
        gtf_check['gene'] = gtf_check['gene'].str.split('\";').str[0]
        gene_instances = len(gtf_check['gene'].unique().tolist())
        if gene_instances == 0:
            raise Exception('No gene or transcript records found in GTF')
        gtf_check = None
    else:
        raise Exception('It appears a properly formatted GTF file was not provided')

    print('Processing', gene_instances, 'gene instances.')
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
                gene_id_original = gtf_remainder.at[end + n, gtf_annotation_column][(gtf_remainder.at[end + n, gtf_annotation_column].find('gene_id \"') + 9):].split('\";')[0]
                gene_id_next = gene_id_original

                amount = -1 # Set scan direction as reverse
                while gene_id_next == gene_id_original:
                    if end + n - 1 < start:
                        amount = 1 # Set scan direction as forward if returned to start gene
                        n = -1 # Reset indexer
                        gtf_remainder = gtf # Point back to original GTF to move past original end set for sub-GTF

                    n += amount # Take another step back until a new gene_id is found (or forward if hit the start)
                    gene_id_next = gtf_remainder.at[end + n, gtf_annotation_column][(gtf_remainder.at[end + n, gtf_annotation_column].find('gene_id \"') + 9):].split('\";')[0]

            if amount == 1: # Correct for scan direction pivot
                n -= 1

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
    ucsc_formatted=False,
    output=True, # True will output all intermediates, not possible if inputting a GTF as pandas dataframe
    threads=None): # Give int for core threshold if desired

    # Import GTF reference file
    if isinstance(gtf, pd.DataFrame) and len(gtf.columns) == 9:
        output = False # Turn off intermediates output
        file_name = None
    elif str(gtf).endswith('.gtf'):
        file_name = str(gtf[:-gtf_rightCoordinate_column]) + '_' # Get rid of GTF extension for now
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
    if longest_transcript == True and ucsc_formatted == False:
        chunks = run_chunks(
            longest_transcripts,
            chunks,
            target_message = 'longest transcripts')

        if output == True:
            file_name = str(file_name) + 'L'
    elif longest_transcript == True and ucsc_formatted == True:
        print("\nWarning: Cannot perform longest transcripts operation on UCSC-formatted GTFs... Skipping...\n")
    else:
        pass

    # Get only protein coding annotated records
    if protein_coding == True:

        if ucsc_formatted == True:
            chunks = run_chunks(
                get_ucsc_proteins,
                chunks,
                target_message = 'protein coding genes')
        else:
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
            _3prime = _3prime,
            ucsc_formatted = ucsc_formatted)
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
