"""
XPRESStools
A toolkit for navigating and analyzing gene expression datasets
alias: xpresstools

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
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
from multiprocessing import cpu_count, Pool
from functools import partial

"""Parse  GTF dataframe for longest transcript per gene record and keep only those transcript records"""
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
                gtf_record = gtf.loc[gene_start_index:gene_stop_index]

                # Find longest transcript for the current gene record
                transcript_lengths = gtf_record.loc[gtf_record[2] == 'transcript'][4] \
                                    - gtf_record.loc[gtf_record[2] == 'transcript'][3]

                transcript_max_index = transcript_lengths.idxmax()

                # Get index for longest transcript
                longest_transcript_id = gtf_record.at[transcript_max_index, 8][(gtf_record.at[transcript_max_index, 8].find('transcript_id \"') + 15):].split('\";')[0]

                # Append longest transcript records to list of records to keep
                long_transcripts.append(gtf_record.loc[gtf_record[8].str.contains(longest_transcript_id)])

    gtf_longest = pd.concat(long_transcripts)
    gtf_longest = gtf_longest.reset_index(drop=True)

    return gtf_longest

"""Only keep GTF records that are annotated as protein_coding"""
def protein_gtf(
    gtf):

    # Take only records that are annotated as 'protein coding'
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]
    gtf_coding = gtf_coding.reset_index(drop=True)

    return gtf_coding

"""Scan first exons recursively by chromosome position and truncate"""
def scan_forward(
    gtf,
    index,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    _5prime,
    _3prime,
    penalty=0):

    # Forward scan for first exon
    n = 0 + penalty
    item = gtf.at[index, 2]
    while item != str(search_string):

        n += 1
        if (index + n) <= (len(gtf.index) - 1): # Make sure next record won't run out of bounds
            item = gtf.at[index + n, 2]
        else:
            return gtf, bad_exons

        if item == str(stop_string): # Make sure didn't run into next transcript
            return gtf, bad_exons

        # Sanity check that we are at an exon record
        if gtf.at[index + n, 2] == str(search_string) \
        and str(annotation) in gtf.at[index + n, 8]:
            # Exon 1 will be first exon
            if gtf.at[index + n, 6] == '+':
                gtf, bad_exons = plus_5prime(
                                    gtf,
                                    index,
                                    n,
                                    _5prime,
                                    _3prime,
                                    bad_exons,
                                    search_string,
                                    stop_string,
                                    annotation,
                                    penalty)
                return gtf, bad_exons

            # Exon 1 will be last exon
            elif gtf.at[index + n, 6] == '-':
                gtf, bad_exons = minus_3prime(
                                    gtf,
                                    index,
                                    n,
                                    _5prime,
                                    _3prime,
                                    bad_exons,
                                    search_string,
                                    stop_string,
                                    annotation,
                                    penalty)
                return gtf, bad_exons

            # Exon 1 does not have a strandedness annotation
            else:
                raise Exception('Unstranded transcript record present')

"""Truncate 5' amount from the first listed positive stranded exon"""
def plus_5prime(
    gtf,
    index,
    counter,
    _5prime,
    _3prime,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[index + counter, 3] + _5prime <= gtf.at[index + counter, 4]:
        gtf.at[index + counter, 3] = gtf.at[index + counter, 3] \
                                    + _5prime
        return gtf, bad_exons

    # Add current exon to list of exons too short and enter recursive loop
    else:
        bad_exons.append(index + counter) # Remove short exon from record
        remainder = _5prime \
                    - abs(gtf.at[index + counter, 4] - gtf.at[index + counter, 3]) # Take what's left over
        return scan_forward( # Recursive scan to next exon until no remainder
                gtf,
                index,
                bad_exons,
                search_string,
                stop_string,
                annotation,
                remainder,
                _3prime,
                penalty + 1)

"""Truncate 3' amount from the first listed minus stranded exon"""
def minus_3prime(
    gtf,
    index,
    counter,
    _5prime,
    _3prime,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[index + counter, 3] + _3prime <= gtf.at[index + counter, 4]:
        gtf.at[index + counter, 3] = gtf.at[index + counter, 3] \
                                        + _3prime
        return gtf, bad_exons

    # Add current exon to list of exons too short and enter recursive loop
    else:
        bad_exons.append(index + counter) # Remove short exon from record
        remainder = _3prime \
                    - abs(gtf.at[index + counter, 4] - gtf.at[index + counter, 3]) # Take what's left over
        return scan_forward( # Recursive scan to next exon until no remainder
                gtf,
                index,
                bad_exons,
                search_string,
                stop_string,
                annotation,
                _5prime,
                remainder,
                penalty + 1)

"""Scan last exons recursively by chromosome position and truncate"""
def scan_backward(
    gtf,
    index,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    _5prime,
    _3prime,
    penalty=0):

    n = 0
    item = ''
    while item != str(stop_string):

        n += 1
        if (index + n) <= (len(gtf.index) - 1): # Make sure next selection not out of bounds
            item = gtf.at[index + n, 2]
        else:
            return gtf, bad_exons

        if item == str(stop_string) \
        or (index + n) == (len(gtf.index) - 1): # Check next selection
            # If next selection is last index and exon, modify coordinates
            y = 0
            item = gtf.at[index + n + y, 2]
            if (index + n) == (len(gtf.index) - 1) \
            and gtf.at[index + n + y, 2] == str(search_string) \
            and str(annotation) in gtf.at[index + n + y, 8]:
                # Exon 1 will be first exon
                if gtf.at[index + n + y, 6] == '+':
                    gtf, bad_exons = plus_3prime(
                                        gtf,
                                        index,
                                        n,
                                        y,
                                        _5prime,
                                        _3prime,
                                        bad_exons,
                                        search_string,
                                        stop_string,
                                        annotation,
                                        penalty)
                    return gtf, bad_exons
                # Exon 1 will be last exon
                elif gtf.at[index + n + y, 6] == '-':
                    gtf, bad_exons = minus_5prime(
                                        gtf,
                                        index,
                                        n,
                                        y,
                                        _5prime,
                                        _3prime,
                                        bad_exons,
                                        search_string,
                                        stop_string,
                                        annotation,
                                        penalty)
                    return gtf, bad_exons
                else:
                    raise Exception('Unstranded transcript record present')

            # Otherwise start tracking back to last exon for transcript
            else:
                y = 0 - penalty
                while item != str(search_string):
                    y -= 1
                    if (index + n + y) > 0: # Check that step back will be valid
                        item = gtf.at[index + n + y, 2]
                    else:
                        return gtf, bad_exons

                    if item == str(stop_string): # Make sure didn't run back into last transcript
                        return gtf, bad_exons

                    # Sanity check that we are at an exon record
                    if gtf.at[index + n + y, 2] == str(search_string) \
                    and str(annotation) in gtf.at[index + n + y, 8]:

                        # Exon last will be last exon
                        if gtf.at[index + n + y, 6] == '+':
                            gtf, bad_exons = plus_3prime(
                                                gtf,
                                                index,
                                                n,
                                                y,
                                                _5prime,
                                                _3prime,
                                                bad_exons,
                                                search_string,
                                                stop_string,
                                                annotation,
                                                penalty)
                            return gtf, bad_exons
                        # Exon last will be first exon
                        elif gtf.at[index + n + y, 6] == '-':
                            gtf, bad_exons = minus_5prime(
                                                gtf,
                                                index,
                                                n,
                                                y,
                                                _5prime,
                                                _3prime,
                                                bad_exons,
                                                search_string,
                                                stop_string,
                                                annotation,
                                                penalty)
                            return gtf, bad_exons
                        else:
                            raise Exception('Unstranded transcript record present')
                    else:
                        pass

"""Truncate 3' amount from the first listed positive stranded exon"""
def plus_3prime(
    gtf,
    index,
    counter,
    inner_counter,
    _5prime,
    _3prime,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[index + counter + inner_counter, 4] - _3prime >= gtf.at[index + counter + inner_counter, 3]:
        gtf.at[index + counter + inner_counter, 4] = gtf.at[index + counter + inner_counter, 4] \
                                                    - _3prime
        return gtf, bad_exons

    # Add current exon to list of exons too short and enter recursive loop
    else:
        bad_exons.append(index + counter + inner_counter) # Remove short exon from record
        remainder = _3prime \
                    - abs(gtf.at[index + counter + inner_counter, 4] - gtf.at[index + counter + inner_counter, 3]) # Take what's left over
        return scan_backward( # Recursive scan to next exon until no remainder
                gtf,
                index + counter + inner_counter,
                bad_exons,
                search_string,
                stop_string,
                annotation,
                _5prime,
                remainder,
                penalty + 1)

"""Truncate 5' amount from the first listed minus stranded exon"""
def minus_5prime(
    gtf,
    index,
    counter,
    inner_counter,
    _5prime,
    _3prime,
    bad_exons,
    search_string,
    stop_string,
    annotation,
    penalty):

    # Edit location and exit the recursive loop
    if gtf.at[index + counter + inner_counter, 4] - _5prime >= gtf.at[index + counter + inner_counter, 3]:
        gtf.at[index + counter + inner_counter, 4] = gtf.at[index + counter + inner_counter, 4] \
                                                    - _5prime
        return gtf, bad_exons

    # Add current exon to list of exons too short and enter recursive loop
    else:
        bad_exons.append(index + counter + inner_counter) # Remove short exon from record
        remainder = _5prime \
                    - abs(gtf.at[index + counter + inner_counter, 4] - gtf.at[index + counter + inner_counter, 3]) # Take what's left over
        return scan_backward( # Recursive scan to next exon until no remainder
                gtf,
                index + counter + inner_counter,
                bad_exons,
                search_string,
                stop_string,
                annotation,
                remainder,
                _3prime,
                penalty + 1)

"""Run MAIN function for GTF truncation"""
def truncate_gtf(
    gtf,
    _5prime=45,
    _3prime=15):

    # Initialize
    gtf_c = gtf.copy() # Make copy in order to edit dataframe

    bad_exons = [] # Make list of indicies with bad exons (too short)

    for index, row in gtf.iterrows():

        # Find records for transcripts
        if row[2] == 'transcript':

            # Recursively scan forward in the transcript to truncate n nucleotides
            gtf_c, bad_exons = scan_forward(gtf_c, index, bad_exons, 'exon', 'transcript', 'exon_number \"', _5prime, _3prime)

            # Forward scan for next transcript, then backtrack to last exon record for transcript
            gtf_c, bad_exons = scan_backward(gtf_c, index, bad_exons, 'exon', 'transcript', 'exon_number \"', _5prime, _3prime)

    # Drop exons that are completely truncated
    print(str(len(bad_exons)) + ' exons records removed from reference chunk for being too short.')

    gtf_c = gtf_c.drop(gtf_c.index[bad_exons])

    return gtf_c

"""Run all GTF-editing functions"""
def edit_gtf(
    gtf, # Dataframe of file path and name to GTF reference
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=True,
    _5prime=45, # If no 5' truncation desired, set to 0
    _3prime=15, # If no 3' truncation desired, set to 0
    output=True, # True will output all intermediates, not possible if inputting a GTF as pandas dataframe
    cpu_threshold=None): # Give int for core threshold if desired

    # Import GTF reference file
    if isinstance(gtf, pd.DataFrame):
        output = False # Turn off intermediates output

        if len(gtf.columns) == 9: # Check some basics on dataframe format
            pass
        else:
            raise Exception('Error: A GTF-formatted dataframe was not provided')

    if str(gtf).endswith('.gtf'):
        file_name = gtf[:-4] # Get rid of GTF extension for now
        gtf = pd.read_csv(str(gtf), sep='\t', header=None, comment='#', low_memory=False)
    else:
        raise Exception('Error: A GTF-formatted file was not provided')

    # Get chunking params
    cores = cpu_count() # Number of CPU cores on your system
    if cpu_threshold == None or cpu_threshold >= cores:
        pass
    else:
        cores = int(cpu_threshold)

    start = 0 # Get first start coordinate for chunk
    batch = round(len(gtf.index) / cores) # Approx. number of samples in a chunk

    # Get chunking indices
    chunks = [] # Initialize chunk storage
    for y in range(cores):

        end = start + batch # Set tentative end

        if end > len(gtf.index) - 1: # If end of dataframe, end there
            end = len(gtf.index) - 1

        else: # To to tentative end and search until next gene record
            gtf_remainder = gtf.iloc[end:]

            gene_id_original = gtf_remainder.at[end, 8][(gtf_remainder.at[end, 8].find('gene_id \"') + 9):].split('\";')[0]
            gene_id_next = gene_id_original
            n = 0

            while gene_id_next == gene_id_original:

                if end + n + 1 > len(gtf.index) - 1:
                    break
                else:
                    n += 1
                    gene_id_next = gtf_remainder.at[end + n, 8][(gtf_remainder.at[end + n, 8].find('gene_id \"') + 9):].split('\";')[0]

        # Parse out current gene record
        chunks.append(gtf.loc[start:end + n])
        start = end + n + 1 # End coordinate for last chunk to start with next

    print('Dataframe split into ' + str(len(chunks)) + ' chunks for parallelization')

    # Run GTF modifications
    # Parse each gene record for longest transcript
    if longest_transcript == True:
        print('Parsing record for longest isoform')
        pool = Pool(cores)
        chunks = pool.map(longest_transcripts, chunks)

        pool.close()
        pool.join()

        if output == True:
            file_name = str(file_name) + '_longestTranscripts'

    # Get only protein coding annotated records
    if protein_coding == True:
        print('Parsing record for protein coding genes')
        pool = Pool(cores)
        chunks = pool.map(protein_gtf, chunks)

        pool.close()
        pool.join()

        if output == True:
            file_name = str(file_name) + '_proteinCoding'

    # Truncate by unique transcript
    # If file has not been parsed for longest transcript per gene, will truncate each isoform
    if truncate_reference == True:
        print('Truncating transcript records')
        pool = Pool(cores)

        func = partial(
                truncate_gtf,
                _5prime = _5prime,
                _3prime = _3prime)
        chunks = pool.map(func, chunks)

        pool.close()
        pool.join()

        if output == True:
            file_name = str(file_name) + '_truncated'
            gtf = pd.concat(chunks)

    # Merge final GTF from chunks and output
    gtf = pd.concat(chunks)

    gtf.to_csv(
        str(file_name) + '.gtf',
        sep = '\t',
        header = None,
        index = False,
        quoting = csv.QUOTE_NONE)
