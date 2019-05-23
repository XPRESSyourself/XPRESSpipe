# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
#__path__, xpresspipe_arguments  =  os.path.split(__file__)
__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'


"""gtfModify functions"""

# Input file will test chunking to make sure GTF structure is conserved
gtf = pd.read_csv(
    str(__path__) + 'other/gtfFlat_test.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)

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


# Test gtfFlatten functions
from xpresspipe.gtfFlatten import flatten_reference

#df = flatten_reference(str(__path__) + 'se_reference/transcripts.gtf')
#df.head()













from xpresspipe.gtfFlatten import create_chromosome_index, create_coordinate_index, flat_list, get_coding_length, make_flatten, flatten_reference
