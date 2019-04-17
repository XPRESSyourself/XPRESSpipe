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
import numpy as np

"""Retrieve meta-coordinate for each read in input BAM matrix"""
def meta_coordinates(
        bam):

    # Map middle point of each read as left-most position plus half of read length
    bam[16] = bam[3] + (bam[9].str.len() / 2).apply(np.floor).astype('int64')

    # Return as array
    mid_coordinates = bam[[2,16]]
    return mid_coordinates.values.tolist()

"""P site ranges for 28-30mers"""
def psite_ranges(
        bam):

    # Keep only optimal footprint size
    bam = bam[(bam[9].str.len() == 28) | (bam[9].str.len() == 29) | (bam[9].str.len() == 30)]

    # Get rightmost coordinates for each read
    bam[16] = bam[3] + bam[9].str.len()

    # Return as array
    phased_coordinates = bam[[2,3,16]]

    if len(phased_coordinates.index) == 0:
        return None
    else:
        return phased_coordinates.values.tolist()


"""Search for plausible records per provided coordinate"""
def get_coordinate_records_meta(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate_start):

    record_array = []
    chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
    for index, record in enumerate(chromosome_array):
        if record[0] <= search_coordinate_start and record[1] >= search_coordinate_start:
            record_array.append(record)

    return record_array

"""Search for plausible records per provided coordinate"""
def get_coordinate_records_period(
        coordinate_index,
        chromosome_index,
        search_chromosome,
        search_coordinate_reverse,
        search_coordinate_forward):

    search_coordinate_reverse += 16
    search_coordinate_forward -= 16

    record_array = []
    chromosome_array = coordinate_index[chromosome_index[search_chromosome]]
    for index, record in enumerate(chromosome_array):
        if record[0] <= search_coordinate_reverse and record[1] >= search_coordinate_reverse \
        or record[0] <= search_coordinate_forward and record[1] >= search_coordinate_forward:
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
    metagene_profile = pd.DataFrame(
        0,
        index = range(101),
        columns = ['metacount'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records_meta(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            position_count = []
            for index, transcript_record in enumerate(record_array):
                position = get_position(
                    record[1],
                    transcript_record[3],
                    transcript_record[2])
                if position != None:
                    position_count.append([position, transcript_record[4]])

            for x in position_count:
                meta_position = int((x[0] / x[1]) * 100)
                count = 1 / len(position_count)
                metagene_profile.at[meta_position, 'metacount'] += count

    return metagene_profile

"""Get periodicity profile for bam file"""
def get_periodicity_profile(
        aligned_reads_index,
        coordinate_index,
        chromosome_index):

    # Initialize profile dataframe for storage
    metagene_profile = pd.DataFrame(
        0,
        index = range(201),
        columns = ['phasing'])

    # Search through each mapped read coordinate
    for index, record in enumerate(aligned_reads_index):
        record_array = get_coordinate_records_period(
                coordinate_index,
                chromosome_index,
                record[0],
                record[1],
                record[2])

        # If a record array is not None, get the exonic position from start for each record for the coordinate
        if record_array:
            position_count = []
            for index, transcript_record in enumerate(record_array):
                # Determine P-site coordinate to search
                if transcript_record[1] == '+':
                    coordinate = record[2] - 16
                else:
                    coordinate = record[1] + 16

                # Get position relative to start of the p-site
                position = get_position(
                    coordinate,
                    transcript_record[3],
                    transcript_record[2])
                if position != None and position < 201:
                    position_count.append(position)

            for x in position_count:
                count = 1 / len(position_count)
                metagene_profile.at[x, 'phasing'] += count

    return metagene_profile
