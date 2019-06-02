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

"""Scan first exons recursively by chromosome position and truncate"""
def scan_forward(
    gtf,
    bad_exons,
    start_index,
    stop_index,
    trim_type,
    _5prime,
    _3prime,
    penalty=0):

    sign = gtf.at[start_index, 6]

    # Start at start_index and parse forward to first exon and trim
    # until stop index (inclusive)
    n = start_index + penalty
    while n != stop_index + 1:

        if gtf.at[n, 2] != trim_type:
            n += 1
            penalty += 1
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
                    penalty)

                return gtf, bad_exons

            # Exon 1 will be last exon
            elif sign == '-':
                gtf, bad_exons = minus_3prime(
                    gtf,
                    bad_exons,
                    n,
                    start_index,
                    stop_index,
                    trim_type,
                    _5prime,
                    _3prime,
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
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[current_position, 3] + _5prime < gtf.at[current_position, 4]:
        gtf.at[current_position, 3] = gtf.at[current_position, 3] \
                                    + _5prime

        return gtf, bad_exons

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:
        bad_exons.append(current_position) # Remove short exon from record
        remainder = _5prime \
                    - abs(gtf.at[current_position, 4] - gtf.at[current_position, 3])

        # Recursive scan to next exon until no remainder
        return scan_forward(
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            remainder,
            _3prime,
            penalty + 1)

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
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[current_position, 3] + _3prime < gtf.at[current_position, 4]:
        gtf.at[current_position, 3] = gtf.at[current_position, 3] \
                                        + _3prime
        return gtf, bad_exons

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:
        bad_exons.append(current_position) # Remove short exon from record
        remainder = _3prime \
                    - abs(gtf.at[current_position, 4] - gtf.at[current_position, 3])

        # Recursive scan to next exon until no remainder
        return scan_forward(
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            _5prime,
            remainder,
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
    penalty=0):

    sign = gtf.at[start_index, 6]

    # Start at start_index and parse forward to first exon and trim
    # until stop index (inclusive)
    n = stop_index - penalty

    while n != start_index - 1:

        if gtf.at[n, 2] != trim_type:
            n -= 1
            penalty -= 1
        else:
            if sign == '+':
                gtf, bad_exons = plus_3prime(
                    gtf,
                    bad_exons,
                    n,
                    start_index,
                    stop_index,
                    trim_type,
                    _5prime,
                    _3prime,
                    penalty)

                return gtf, bad_exons

            # Exon 1 will be last exon
            elif sign == '-':
                gtf, bad_exons = minus_5prime(
                    gtf,
                    bad_exons,
                    n,
                    start_index,
                    stop_index,
                    trim_type,
                    _5prime,
                    _3prime,
                    penalty)

                return gtf, bad_exons

            # Exon 1 does not have a strandedness annotation
            else:
                raise Exception('Unstranded transcript record present')

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
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[current_position, 4] - _3prime > gtf.at[current_position, 3]:
        gtf.at[current_position, 4] = gtf.at[current_position, 4] \
                                                    - _3prime
        return gtf, bad_exons

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:
        bad_exons.append(current_position)
        remainder = _3prime \
                    - abs(gtf.at[current_position, 4] - gtf.at[current_position, 3])

        return scan_backward( # Recursive scan to next exon until no remainder
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            _5prime,
            remainder,
            penalty + 1)

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
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[current_position, 4] - _5prime > gtf.at[current_position, 3]:
        gtf.at[current_position, 4] = gtf.at[current_position, 4] \
                                                    - _5prime
        return gtf, bad_exons

    # Add current exon to list of exons too short, take remainder, and enter recursive loop
    else:
        bad_exons.append(current_position) # Remove short exon from record
        remainder = _5prime \
                    - abs(gtf.at[current_position, 4] - gtf.at[current_position, 3])

        # Recursive scan to next exon until no remainder
        return scan_backward(
            gtf,
            bad_exons,
            start_index,
            stop_index,
            trim_type,
            remainder,
            _3prime,
            penalty + 1)

"""Run MAIN function for GTF truncation"""
def truncate_gtf(
    gtf,
    _5prime=45,
    _3prime=15,
    trim_type='CDS'):

    """STEP 1"""
    # Initialize variables for step 1
    gtf_c = gtf.copy()
    bad_transcript = []
    bad_exons = []
    limit = _5prime + _3prime
    parse_type = 'transcript_id \"'

    """Step 1"""
    # Remove transcripts with exon space smaller than truncation sum
    for index, row in gtf.iterrows():

        # Start at a transcript
        if row[2] == 'transcript':

            # Get transcript ID
            transcript_id = gtf.at[index, 8][(gtf.at[index, 8].find(parse_type) + 15):].split('\";')[0]

            # Find range for this transcript and parse out
            n = 0
            parse_id = transcript_id
            while transcript_id == parse_id:
                n += 1
                if index + n > gtf.index[-1]: # Check if we reached the end of the GTF
                    break
                else:
                    parse_id = gtf.at[index + n, 8][(gtf.at[index + n, 8].find(parse_type) + 15):].split('\";')[0]

            gtf_parse = gtf.loc[index:index + n - 1]

            # Create exon length array for each exon labeled record for the transcript
            exon_lengths = gtf_parse.loc[gtf_parse[2] == trim_type][4] \
                                - gtf_parse.loc[gtf_parse[2] == trim_type][3]

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
                    _3prime)

                # Forward scan for next transcript, then backtrack to last exon record for transcript
                gtf_c, bad_exons = scan_backward(
                    gtf_c,
                    bad_exons,
                    index,
                    index + n - 1,
                    trim_type,
                    _5prime,
                    _3prime)

    """STEP 3"""
    # Do final check that no column 3 value is greater or equal to column 4
    # If the case, append to list and drop
    remaining_bad = []
    for index, row in gtf_c.iterrows():
        if row[3] >= row[4]:
            remaining_bad.append(index)

    """STEP 4"""
    # Remove bad transcripts and reindex GTF
    remove_indices = list(set(bad_exons + bad_transcript + remaining_bad))

    print(str(len(remove_indices)) + ' records removed from reference chunk for being too short.')
    gtf_c = gtf_c.drop(gtf_c.index[remove_indices])
    gtf_c = gtf_c.reset_index(drop=True)

    # Garbage management
    #gtf = None
    #gc.collect()

    return gtf_c
