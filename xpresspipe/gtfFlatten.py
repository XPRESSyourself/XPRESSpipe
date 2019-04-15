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
import pandas as pd

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

"""Search for plausible records per provided coordinate"""
def get_coordinate_records(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate):

    record_array = []
    chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
    for index, record in enumerate(chromosome_array):
        if record[0] <= search_coordinate and record[1] >= search_coordinate:
            record_array.append(record)

    return record_array

"""Get relative position of coordinate to the start of the transcript"""
def get_position(
        position,
        coordinates,
        strand):

    # Can order operations same since reference puts - strand records in reverse already
    location = 0
    last_coordinate = 0

    for y in coordinates:

        # Exit if mapped to intron
        if strand == '+':
            next_coordinate = min(y)
        else: # '-'
            next_coordinate =  max(y)

        if position > last_coordinate and position < next_coordinate:
            return None

        # Map to exon position
        if position >= min(y) and position <= max(y):
            if strand == '+':
                location += abs(position - min(y))
            else: # '-'
                location += abs(max(y)- position)

            return location

        else:
            location += abs(y[1] - y[0])
            if strand == '+':
                next_coordinate = max(y)
            else: # '-'
                next_coordinate =  min(y)

"""Get meta profile for bam file"""
def get_meta_profile(
        aligned_reads_index,
        coordinate_index,
        chromosome_index):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(0, index=range(101), columns=['metacount'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            for index, transcript_record in enumerate(record_array):
                position = get_position(
                    record[1],
                    transcript_record[3],
                    transcript_record[2])
                if position != None:
                    meta_position = int((position / transcript_record[4]) * 100)
                    count = 1 / len(record_array)
                    metagene_profile.at[meta_position, 'metacount'] += count

    return metagene_profile

"""Get periodicity profile for bam file"""
def get_periodicity_profile(
        aligned_reads_index,
        coordinate_index,
        chromosome_index):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(0, index=range(201), columns=['metacount'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            for index, transcript_record in enumerate(record_array):
                position = get_position(
                    record[1],
                    transcript_record[3],
                    transcript_record[2])
                if position != None and position < 201:
                    count = 1 / len(record_array)
                    metagene_profile.at[position, 'metacount'] += count

    return metagene_profile

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
        transcript_length += abs(int(y[1]) - int(y[0]))

    return transcript_length

"""Make a flat reference from parsed GTF file"""
def make_flatten(
        gtf):

    records = []

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

                    if gtf.at[index + n, 2] == 'exon': # Append coordinate paires for each exon of the transcript
                        coordinates.append([gtf.at[index + n, 3], gtf.at[index + n, 4]])

            # Get start and end positions for transcript/gene
            # Assumes start and stop will be start and end of a protein coding transcript
            start = min(flat_list(coordinates))
            end = max(flat_list(coordinates))

            # Push information to record
            records.append([gene, strand, chromosome, start, end, coordinates])

    # Push flattened reference into pandas dataframe
    headers = ['gene', 'strand', 'chromosome', 'start', 'end', 'coordinates']
    reference = pd.DataFrame(records, columns=headers)
    reference.chromosome = reference.chromosome.astype(str)

    # Get length of each transcripts exon space
    reference['length'] = reference['coordinates'].apply(get_coding_length)

    return reference

"""Read in coordinate reference from GTF"""
def flatten_reference(
        gtf,
        threads=None):

    # Import GTF reference file
    if isinstance(gtf, pd.DataFrame) and len(gtf.columns) == 9:
        pass
    elif str(gtf).endswith('.gtf'):
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
        gtf = gtf.reset_index(drop=True)
        gtf = make_flatten(gtf)

        return gtf

    else:
        return
