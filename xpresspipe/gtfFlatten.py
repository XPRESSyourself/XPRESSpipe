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
import pandas as pd
import gc

"""IMPORT INTERNAL DEPENDENCIES"""
from .gtfModify import get_chunks, run_chunks, longest_transcripts, protein_gtf

"""Get index number for each chromosome name"""
def create_chromosome_index(
    gtf_flat):

    chromosome_index = {}
    y = 0
    for x in list(gtf_flat['chromosome'].unique()):
        chromosome_index[x] = y
        y += 1

    return chromosome_index

"""Create a chromosome indexed array with positions and key"""
def create_coordinate_index(
    gtf_flat):

    coordinate_index = []
    for x in list(gtf_flat['chromosome'].unique()):
        coordinate_index.append(gtf_flat[['start','end','strand','coordinates','length']][gtf_flat['chromosome'] == str(x)].values.tolist())

    return coordinate_index

"""Flatten nested list"""
def flat_list(
    nested_list):

    flat_list = []
    for x in nested_list:
        for y in x:
            flat_list.append(y)

    return flat_list

"""Get coding space length for each transcript"""
def get_coding_length(
    coordinates):

    transcript_length = 0
    for y in coordinates:
        transcript_length += abs(int(y[1]) - int(y[0])) + 1

    return transcript_length

"""Make a flat reference from parsed GTF file"""
def make_flatten(
    gtf,
    record_type='exon'):

    print('Flattening reference file...')

    records = []
    warning = 0
    for index, row in gtf.iterrows():

        if row[2] == 'transcript':

            gene = gtf.at[index, 8][(gtf.at[index, 8].find('gene_id \"') + 9):].split('\";')[0]
            strand = row[6]
            chromosome = row[0]
            coordinates = []

            # Scan for exons for the given
            n = 0
            item = gene
            while item == gene:
                n += 1

                if index + n > len(gtf.index) - 1 or gtf.at[index + n, 2] == 'transcript':
                    break
                else:
                    item = gtf.at[index + n, 8][(gtf.at[index + n, 8].find('gene_id \"') + 9):].split('\";')[0]

                    if gtf.at[index + n, 2] == record_type.lower(): # Append coordinate paires for each exon of the transcript
                        coordinates.append([gtf.at[index + n, 3], gtf.at[index + n, 4]])

            # Get start and end positions for transcript/gene
            # Assumes start and stop will be start and end of a protein coding transcript
            try:
                start = min(flat_list(coordinates))
                end = max(flat_list(coordinates))

                # Push information to record
                records.append([gene, strand, chromosome, start, end, coordinates])
            except:
                warning = 1
                print('Warning: No ' + str(record_type.lower()) +  ' records found for gene record ' + str(gene))

    # Push flattened reference into pandas dataframe
    gtf = None
    headers = ['gene', 'strand', 'chromosome', 'start', 'end', 'coordinates']
    reference = pd.DataFrame(records, columns=headers)
    records = None
    gc.collect()
    reference.chromosome = reference.chromosome.astype(str)

    # Get length of each transcripts exon space
    reference['length'] = reference['coordinates'].apply(get_coding_length)

    # Print message about warning
    if warning == 1:
        print('If no gene records were found for a given gene, it is likely in the GTF file parsing where no suitable record was found ')

    return reference

"""
Read in coordinate reference from GTF
Requires an unmodified GTF file
"""
def flatten_reference(
    gtf_file,
    record_type='exon',
    threads=None):

    # Import GTF reference file
    if str(gtf_file).endswith('.gtf'):
        gtf = pd.read_csv(
            str(gtf_file),
            sep = '\t',
            header = None,
            comment = '#',
            low_memory = False)
    elif isinstance(gtf_file, pd.DataFrame):
        gtf = gtf_file

    else:
        raise Exception('Error: A GTF-formatted file or dataframe was not provided')

    gtf_file = None
    del gtf_file
    gc.collect()

    # Get chunks
    chunks = get_chunks(
        gtf,
        threads = threads)

    # Get longest transcripts
    chunks = run_chunks(
        longest_transcripts,
        chunks,
        target_message = 'longest transcripts')

    # Get only protein coding annotated records
    chunks = run_chunks(
        protein_gtf,
        chunks,
        target_message = 'protein coding genes')

    # Rejoin chunks into single GTF and flatten
    if len(chunks) > 0:
        gtf = pd.concat(chunks)
        chunks = None
        del chunks
        gtf = gtf.reset_index(drop=True)
        gtf = make_flatten(gtf, record_type)

        return gtf

    else:
        return
