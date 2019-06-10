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
import re
import sys

"""Get overrepresented sequences from a given FastQC zip file"""
def get_overrep_seqs(directory):

    file = str(directory) + '/fastqc_data.txt'
    if os.path.isfile(file) == False:
        raise Exception('No fastqc data file found in fastqc folder: ' + str(directory))


    # If execution makes it to here, should have a valid data filename
    with folder.open(file) as datafile:
        overrep_seqs = []
        store_flag = False
        for bytes_string in datafile:
            line = str(bytes_string, 'utf-8').strip()
            if line == ">>Overrepresented sequences	fail":
                # This means time to start storing lines
                store_flag = True
            elif store_flag and line == ">>END_MODULE":
                break
            elif store_flag:
                entry = line.split()[:2]
                try:
                    entry[1] = int(entry[1])
                    overrep_seqs.append(entry)
                except Exception as e:
                    pass

    if len(overrep_seqs) <= 0:
        print ("No fastqc data file found in fastqc folder: " + file)
        sys.exit(1)
    return overrep_seqs

"""Check for a match of the sequence on right side of new with anywhere in old"""
def getMatchRight(seq, freq, combined_seqs, min_ovlp):
    # Take the rightmost chars from the new seq
    match_chars = seq[-min_ovlp:]
    found = False
    to_remove = -1
    to_add = []
    for i in range(len(combined_seqs)):
        old_seq = combined_seqs[i][0]
        idx = old_seq.find(match_chars)
        if idx != -1:
            poss_match_len = idx+min_ovlp
            if poss_match_len > len(seq):
                poss_match_len = len(seq)
            if old_seq[:idx+min_ovlp] == seq[-poss_match_len:]:
                # It's a match so remove old entry and create a new one with the combined seq+freq
                to_add = [seq[:-poss_match_len] + old_seq, freq+combined_seqs[i][1]]
                to_remove = i
                found = True
    if found:
        combined_seqs.pop(to_remove)
    return combined_seqs,to_add,found

"""Check for a match of the sequence on left side of new with anywhere in old"""
def getMatchLeft(seq, freq, combined_seqs, min_ovlp):
    # Take the leftmost chars from the new seq
    match_chars = seq[:min_ovlp]
    found = False
    to_remove = -1
    to_add = []
    for i in range(len(combined_seqs)):
        old_seq = combined_seqs[i][0]
        idx = old_seq.rfind(match_chars)
        if idx != -1:
            poss_match_end = 0
            poss_match_len = len(old_seq)-idx
            if poss_match_len > len(seq):
                poss_match_len = len(seq)
            if old_seq[idx:poss_match_len+idx] == seq[:poss_match_len]:
                to_add = [old_seq+seq[poss_match_len:], freq+combined_seqs[i][1]]
                to_remove = i
                found = True
    if found:
        combined_seqs.pop(to_remove)
    return combined_seqs,to_add,found

"""Add entry to list of overrepresented sequences and count"""
def addEntry(new_entry, combined_entries, min_ovlp):
    new_seq = new_entry[0]
    new_freq = new_entry[1]
    found = False
    to_add = []
    for i in range(len(combined_entries)):
        old_seq = combined_entries[i][0]
        # The new seq is a subset of the old one, so combine freqs
        if new_seq in old_seq:
            combined_entries[i][1] += new_freq
            found = True
        # The old seq is a subset of the new one, so combine freqs
        elif old_seq in new_seq:
            combined_entries[i][1] += new_freq
            combined_entries[i][0] = new_seq
            found = True
        else:
            combined_entries,to_add,found = getMatchRight(new_seq, new_freq, combined_entries, min_ovlp)
            if not found:
                combined_entries,to_add,found = getMatchLeft(new_seq, new_freq, combined_entries, min_ovlp)
        if found:
            break
    if found and len(to_add) > 0:
        addEntry(to_add, combined_entries,min_ovlp)
    elif not found:
        combined_entries.append(new_entry)
    return combined_entries

def combineSeqs(entries, min_overlap):
    combined_entries = []
    for entry in entries:
        combined_entries = addEntry(entry, combined_entries, min_overlap)
    return combined_entries

def countFreqs(combined_entries):
    count = 0
    for entry in combined_entries:
        count += entry[1]
    return count

"""Determine consensus overrepresented sequences between files from list"""
def rrnaProbe(
    directory_list,
    min_overlap):

    footprint_seqs = []

    for directory in directory_list:
        footprint_seqs += get_overrep_seqs(directory)

    combined_fp = []

    # Put into uppercase format together, prep for combining (ribosome footprint sequences)
    count = 0
    for seq,count in footprint_seqs:
        seq = seq.strip().upper()
        combined_fp.append([seq,count])
    combined_fp = combineSeqs(combined_fp, min_overlap)

    combined_fp.sort(key=lambda x: x[1], reverse=True)
    results_list = []
    for entry in combined_fp:
        results_list.append("\t".join([str(x) for x in entry]))

    header = "##This file contains the most highly represented sequences in a given ribosome footprint RNA-seq experiment\n" +\
            "##The sequences included are combined from FASTQC output so that if multiple sequences reported there are exact-match substring of one another, they will be combined and their counts summed\n"+\
            "##These sequences should be checked using a tool such as BLAST to ensure that the most represented sequences are not Illumina artifacts.\n" +\
            "#SEQ\tCOUNT\n"

    return header + "\n".join(results_list)
