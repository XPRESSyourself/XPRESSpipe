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

"""GLOBALS"""
gtf_chr_column = 1
gtf_type_column = 2
gtf_leftCoordinate_column = 3
gtf_rightCoordinate_column = 4
gtf_sign_column = 6
gtf_annotation_column = 8
parse_type = 'transcript_id \"'

"""Scan first exons recursively by chromosome position and truncate"""
def scan_forward(
    gtf,
    bad_exons,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    #recursive_counter,
    penalty=0):

    sign = gtf.at[start_index, gtf_sign_column]

    # Start at start_index and parse forward to first exon and trim
    # until stop index (inclusive)
    n = start_index + penalty

    if n > stop_index:
        return

    while n != stop_index + 1:

        if gtf.at[n, gtf_type_column] != trim_type:
            n += 1
            penalty += 1 # Make sure penalty keeps up with current position steps
        else:
            if sign == '+':
                gtf, bad_exons = plus_5prime(
                    gtf,
                    bad_exons,
                    n,
                    start_index,
                    stop_index,
                    trim_type,
                    _5prime,
                    _3prime,
                    ucsc_formatted,
                    penalty)

                return gtf, bad_exons

            # Exon 1 will be last exon
            elif sign == '-':
                if ucsc_formatted == True:
                    gtf, bad_exons = plus_5prime(
                        gtf,
                        bad_exons,
                        n,
                        start_index,
                        stop_index,
                        trim_type,
                        _3prime, # swap order for UCSC to make use of plus_5prime with the 3prime amount
                        _5prime,
                        ucsc_formatted,
                        penalty)
                else:
                    gtf, bad_exons = minus_5prime(
                        gtf,
                        bad_exons,
                        n,
                        start_index,
                        stop_index,
                        trim_type,
                        _5prime,
                        _3prime,
                        ucsc_formatted,
                        penalty)

                return gtf, bad_exons

            # Exon 1 does not have a strandedness annotation
            else:
                raise Exception('Unstranded transcript record present')

"""Truncate 5' amount from the first listed positive stranded exon"""
def plus_5prime(
    gtf,
    bad_exons,
    current_position,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    penalty):
    #recursive_counter):

    # Edit location and exit the recursive loop
    if _5prime \
    - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
    - gtf.at[current_position, gtf_leftCoordinate_column]
    + 1) == 0:

        bad_exons.append(current_position) # Remove short exon from record
        return gtf, bad_exons #, recursive_counter

    elif gtf.at[current_position, gtf_leftCoordinate_column] + _5prime \
    <= gtf.at[current_position, gtf_rightCoordinate_column]:

        gtf.at[current_position, gtf_leftCoordinate_column] \
            = gtf.at[current_position, gtf_leftCoordinate_column] \
            + _5prime

        return gtf, bad_exons #, recursive_counter

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:

        bad_exons.append(current_position) # Remove short exon from record
        remainder = _5prime \
            - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
            - gtf.at[current_position, gtf_leftCoordinate_column]
            + 1)

        #recursive_counter += 1

        # Recursive scan to next exon until no remainder
        return scan_forward(
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            remainder,
            _3prime,
            ucsc_formatted,
            #recursive_counter,
            penalty + 1) # Penalty takes current position plus one for next round

"""Truncate 5' amount from the first listed minus stranded exon"""
def minus_5prime(
    gtf,
    bad_exons,
    current_position,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    penalty):
    #recursive_counter):

    # Edit location and exit the recursive loop
    if _5prime \
    - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
    - gtf.at[current_position, gtf_leftCoordinate_column] \
    + 1) == 0:
        bad_exons.append(current_position) # Remove short exon from record
        return gtf, bad_exons #, recursive_counter

    elif gtf.at[current_position, gtf_rightCoordinate_column] - _5prime \
    >= gtf.at[current_position, gtf_leftCoordinate_column]:

        gtf.at[current_position, gtf_rightCoordinate_column] \
            = gtf.at[current_position, gtf_rightCoordinate_column] \
            - _5prime
        return gtf, bad_exons #, recursive_counter

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:

        bad_exons.append(current_position) # Remove short exon from record
        remainder = _5prime \
            - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
            - gtf.at[current_position, gtf_leftCoordinate_column] \
            + 1)

        #recursive_counter += 1

        if ucsc_formatted == True:
            return scan_backward(
                gtf,
                bad_exons,
                start_index,
                stop_index,
                trim_type,
                remainder,
                _3prime,
                ucsc_formatted,
                #recursive_counter,
                penalty + 1)

        else:
            # Recursive scan to next exon until no remainder
            return scan_forward(
                gtf,
                bad_exons,
                start_index,
                stop_index,
                trim_type,
                remainder,
                _3prime,
                ucsc_formatted,
                #recursive_counter,
                penalty + 1)

"""Scan first exons recursively by chromosome position and truncate"""
def scan_backward(
    gtf,
    bad_exons,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    #recursive_counter,
    penalty=0):

    sign = gtf.at[start_index, gtf_sign_column]

    # Start at start_index and parse forward to first exon and trim
    # until stop index (inclusive)
    n = stop_index - penalty

    while n != start_index - 1:
        # Built in control to prevent indexing out of bounds, not relevant for small tests, but required when working with larger GTFs
        if n > gtf.index[-1] or n < start_index:
            return gtf, bad_exons # Says we are out of bounds and return the work that its done up until now

        if gtf.at[n, gtf_type_column] != trim_type:
            n -= 1
            penalty += 1 # Make sure penalty keeps up with current position steps

        else:
            if sign == '+':
                gtf, bad_exons = plus_3prime( # add recursive_counter to outputs
                    gtf,
                    bad_exons,
                    n,
                    start_index,
                    stop_index,
                    trim_type,
                    _5prime,
                    _3prime,
                    ucsc_formatted,
                    penalty)
                    #recursive_counter)

                return gtf, bad_exons #, recursive_counter

            # Exon 1 will be last exon
            elif sign == '-':
                if ucsc_formatted == True:
                    gtf, bad_exons = minus_5prime( # add recursive_counter to outputs
                        gtf,
                        bad_exons,
                        n,
                        start_index,
                        stop_index,
                        trim_type,
                        _5prime,
                        _3prime,
                        ucsc_formatted,
                        penalty)
                        #recursive_counter)
                else:
                    gtf, bad_exons = minus_3prime( # add recursive_counter to outputs
                        gtf,
                        bad_exons,
                        n,
                        start_index,
                        stop_index,
                        trim_type,
                        _5prime,
                        _3prime,
                        ucsc_formatted,
                        penalty)
                        #recursive_counter)

                return gtf, bad_exons #, recursive_counter

            # Exon 1 does not have a strandedness annotation
            else:
                raise Exception('Unstranded transcript record present')
    else:
        return gtf, bad_exons #, recursive_counter # Provide exit mechanism if start hit outside while loop

"""Truncate 3' amount from the first listed positive stranded exon"""
def plus_3prime(
    gtf,
    bad_exons,
    current_position,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    penalty):
    #recursive_counter):

    # Edit location and exit the recursive loop
    if _3prime \
    - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
    - gtf.at[current_position, gtf_leftCoordinate_column] \
    + 1) == 0:

        bad_exons.append(current_position)
        return gtf, bad_exons #, recursive_counter

    elif gtf.at[current_position, gtf_rightCoordinate_column] - _3prime \
    >= gtf.at[current_position, gtf_leftCoordinate_column]:

        gtf.at[current_position, gtf_rightCoordinate_column] \
            = gtf.at[current_position, gtf_rightCoordinate_column] \
            - _3prime

        return gtf, bad_exons #, recursive_counter

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:

        bad_exons.append(current_position)
        remainder = _3prime \
            - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
            - gtf.at[current_position, gtf_leftCoordinate_column] \
            + 1)

        return scan_backward( # Recursive scan to next exon until no remainder
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            _5prime,
            remainder,
            ucsc_formatted,
            #recursive_counter,
            penalty + 1) # Penalty takes current position plus one for next round

"""Truncate 3' amount from the first listed minus stranded exon"""
def minus_3prime(
    gtf,
    bad_exons,
    current_position,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    ucsc_formatted,
    penalty):
    #recursive_counter):

    # Edit location and exit the recursive loop
    if _3prime \
    - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
    - gtf.at[current_position, gtf_leftCoordinate_column] \
    + 1) == 0:

        bad_exons.append(current_position) # Remove short exon from record
        return gtf, bad_exons #, recursive_counter

    elif gtf.at[current_position, gtf_leftCoordinate_column] + _3prime \
    <= gtf.at[current_position, gtf_rightCoordinate_column]:

        gtf.at[current_position, gtf_leftCoordinate_column] \
            = gtf.at[current_position, gtf_leftCoordinate_column] \
            + _3prime

        return gtf, bad_exons #, recursive_counter

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:

        bad_exons.append(current_position) # Remove short exon from record
        remainder = _3prime \
            - abs(gtf.at[current_position, gtf_rightCoordinate_column] \
            - gtf.at[current_position, gtf_leftCoordinate_column] \
            + 1)

        # Recursive scan to next exon until no remainder
        if ucsc_formatted == True:
            return scan_forward(
                gtf,
                bad_exons,
                start_index,
                stop_index,
                trim_type,
                _5prime,
                remainder,
                ucsc_formatted,
                #recursive_counter,
                penalty + 1) # Penalty takes current position plus one for next round
        else:
            return scan_backward(
                gtf,
                bad_exons,
                start_index,
                stop_index,
                trim_type,
                _5prime,
                remainder,
                ucsc_formatted,
                #recursive_counter,
                penalty + 1) # Penalty takes current position plus one for next round

"""Run MAIN function for GTF truncation"""
def truncate_gtf(
    gtf,
    _5prime=45,
    _3prime=15,
    ucsc_formatted=False,
    trim_type='CDS'):

    """STEP 1"""
    # Initialize variables for step 1
    gtf_c = gtf.copy()
    bad_transcript = []
    bad_exons = []
    limit = _5prime + _3prime
    parse_id = ''

    """Step 1"""
    # Remove transcripts with exon space smaller than truncation sum
    for index, row in gtf.iterrows():

        # Get transcript ID
        if parse_type in row[gtf_annotation_column]:
            transcript_id = gtf.at[index, gtf_annotation_column][(gtf.at[index, gtf_annotation_column].find(parse_type) + 15):].split('\";')[0]

            if transcript_id != parse_id:
                parse_id = transcript_id
                # Find range for this transcript and parse out
                n = 0
                while transcript_id == parse_id:
                    n += 1

                    # Check if we reached the end of the GTF
                    if index + n > gtf.index[-1]:
                        break
                    # Handles a GTF where gene record present -- without this
                    # it will try to parse transcript_id incorrectly from row
                    # and screw things up
                    elif parse_type in gtf.at[index + n, gtf_annotation_column]:
                        parse_id = gtf.at[index + n, gtf_annotation_column][(gtf.at[index + n, gtf_annotation_column].find(parse_type) + 15):].split('\";')[0]
                    else:
                        break

                gtf_parse = gtf.loc[index:index + n - 1]
                parse_id = transcript_id # reset for next round check

                # If target type not in transcript record, skip
                # Without, will remove all exons and transcript IDs
                if trim_type not in gtf_parse[gtf_type_column].tolist():
                    continue

                # Create exon length array for each exon labeled record for the transcript
                exon_lengths = gtf_parse.loc[gtf_parse[gtf_type_column] == trim_type][gtf_rightCoordinate_column] \
                    - gtf_parse.loc[gtf_parse[gtf_type_column] == trim_type][gtf_leftCoordinate_column] \
                    + 1

                # Check exons for this record and check
                if exon_lengths.sum() <= limit:

                    # If failing, add transcript_id to bad list
                    for x in range(index, index + n):
                        bad_transcript.append(x)

                else:
                    """STEP 2"""
                    # Run truncator
                    # Recursively scan forward in the transcript to truncate n nucleotides
                    gtf_c, bad_exons = scan_forward(
                        gtf_c,
                        bad_exons,
                        index,
                        index + n - 1,
                        trim_type,
                        _5prime,
                        _3prime,
                        ucsc_formatted)

                    # Forward scan for next transcript, then backtrack to last exon record for transcript
                    gtf_c, bad_exons = scan_backward(
                        gtf_c,
                        bad_exons,
                        index,
                        index + n - 1,
                        trim_type,
                        _5prime,
                        _3prime,
                        ucsc_formatted)

    """STEP 3"""
    # Do final check that no column 3 value is greater or equal to column 4
    # If the case, append to list and drop
    remaining_bad = []
    for index, row in gtf_c.iterrows():
        if row[gtf_leftCoordinate_column] > row[gtf_rightCoordinate_column] \
        and row[gtf_type_column] == trim_type:
            remaining_bad.append(index)

    """STEP 4"""
    # Remove bad transcripts and reindex GTF
    remove_indices = list(set(bad_exons + bad_transcript + remaining_bad))

    print(str(len(remove_indices)) + ' records removed from reference chunk for being too short.')
    #print(str(recursive_counter) + ' transcripts requiring recursive truncation.')
    gtf_c = gtf_c.drop(remove_indices)
    gtf_c = gtf_c.reset_index(drop=True)

    return gtf_c
