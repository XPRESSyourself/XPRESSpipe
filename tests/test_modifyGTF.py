# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
__path__  =  os.path.dirname(os.path.realpath(__file__))
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'

"""gtfModify functions"""

# Input file will test chunking to make sure GTF structure is conserved
gtf = pd.read_csv(
    str(__path__) + 'other/gtf_test.gtf',
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

gtf = pd.read_csv(
    str(__path__) + 'other/gtf_test_large.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)

### Functional tests
from xpresspipe.gtfModify import longest_transcripts, protein_gtf, edit_gtf
gtf_long = longest_transcripts(gtf)
gtf_protein = protein_gtf(gtf_long)
gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=True,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit.shape == (14, 9), 'Something went wrong during the functional test running through all GTF modifications'

print(gtf_edit)
print('----------')
### Unit tests
# Get longest transcript records only for each gene
gtf_long = longest_transcripts(gtf)

# Get protein coding only records
gtf_protein = protein_gtf(gtf_long)

# Run whole gambit together
gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)

gtf_edit = edit_gtf(
    gtf,
    longest_transcript=False,
    protein_coding=True,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)

gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=False,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)

gtf_edit_truncated = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=False,
    truncate_reference=True,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
