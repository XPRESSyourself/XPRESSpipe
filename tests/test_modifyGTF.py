# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
__path__  =  os.path.dirname(os.path.realpath(__file__))
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'

print(__path__)
"""gtfModify functions"""

# Input file will test chunking to make sure GTF structure is conserved
gtf = pd.read_csv(
    str(__path__) + 'other/gtfFlat_test.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)

# Test chunking
from xpresspipe.gtfModify import get_chunks

chunks = get_chunks(gtf, threads=0)
assert len(chunks) == 1, 'get_chunks() failed to catch too few threads'

chunks = get_chunks(gtf, threads=1)
assert len(chunks) == 1, 'get_chunks() failed to limit to 1 thread'

chunks = get_chunks(gtf, threads=2)
assert len(chunks) == 2, 'get_chunks() failed with 2 thread option for test GTF'

chunks = get_chunks(gtf, threads=3)
assert len(chunks) == 2, 'get_chunks() failed to catch too many threads'

chunks = get_chunks(gtf)
assert len(chunks) == 2, 'get_chunks() failed with no thread option provided for test GTF'
assert chunks[0].shape == (31,9) and chunks[1].shape == (33,9), 'get_chunks() failed output proper size chunks'

# Get longest transcript records only for each gene
from xpresspipe.gtfModify import longest_transcripts
gtf_long_truth = [
    [19, 'havana', 'transcript', 60951, 70976],
    [19, 'havana', 'exon', 70928, 70976],
    [19, 'havana', 'exon', 66346, 66499],
    [19, 'havana', 'exon', 60951, 61894],
    [19, 'havana', 'transcript', 107104, 117102],
    [19, 'havana', 'exon', 107104, 107157],
    [19, 'havana', 'exon', 107473, 107555],
    [19, 'havana', 'exon', 110625, 110681],
    [19, 'havana', 'exon', 116507, 117102],]
gtf_long_truth = pd.DataFrame(gtf_long_truth)
gtf_long = longest_transcripts(gtf)
assert gtf_long.iloc[:,:5].equals(gtf_long_truth), 'longest_transcripts() failed'

# Get protein coding only records
from xpresspipe.gtfModify import protein_gtf
gtf_protein_truth = [
    [19, 'havana', 'transcript', 107104, 117102],
    [19, 'havana', 'exon', 107104, 107157],
    [19, 'havana', 'exon', 107473, 107555],
    [19, 'havana', 'exon', 110625, 110681],
    [19, 'havana', 'exon', 116507, 117102],]
gtf_protein_truth = pd.DataFrame(gtf_protein_truth)
gtf_protein = protein_gtf(gtf_long)
assert gtf_protein.iloc[:,:5].equals(gtf_protein_truth), 'protein_gtf() failed'

# Run whole gambit together
from xpresspipe.gtfModify import edit_gtf
gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit.iloc[:,:5].equals(gtf_protein_truth), 'edit_gtf() failed'

gtf_edit = edit_gtf(
    gtf,
    longest_transcript=False,
    protein_coding=True,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit.shape == (33,9), 'edit_gtf() failed during protein_coding only output'

gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=False,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit.shape == (9,9), 'edit_gtf() failed during longest_transcript only output'

gtf_trunc_truth = [
    [19, 'havana', 'transcript', 60951, 70976],
    [19, 'havana', 'exon', 70943, 70976],
    [19, 'havana', 'exon', 66346, 66499],
    [19, 'havana', 'exon', 60951, 61849],
    [19, 'havana', 'transcript', 107104, 117102],
    [19, 'havana', 'exon', 107149, 107157],
    [19, 'havana', 'exon', 107473, 107555],
    [19, 'havana', 'exon', 110625, 110681],
    [19, 'havana', 'exon', 116507, 117087],]
gtf_trunc_truth = pd.DataFrame(gtf_trunc_truth)
gtf_edit_truncated = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=False,
    truncate_reference=True,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit_truncated.iloc[:,:5].equals(gtf_trunc_truth), 'edit_gtf() failed during truncation'
