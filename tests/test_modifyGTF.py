# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
from multiprocessing import cpu_count
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

"""gtfModify functions"""

# Input file will test chunking to make sure GTF structure is conserved
gtf = pd.read_csv(
    str(__path__) + 'references/other/gtf_test.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)
print('GTF Shape:',gtf.shape)

# Test chunking
from xpresspipe.gtfModify import get_chunks

chunks = get_chunks(gtf, threads=0)
assert len(chunks) == 1, 'get_chunks() failed to catch too few threads'

chunks = get_chunks(gtf, threads=1)
assert len(chunks) == 1, 'get_chunks() failed to limit to 1 thread'

chunks = get_chunks(gtf, threads=2)
assert len(chunks) == 2, 'get_chunks() failed with 2 thread option for test GTF'

chunks = get_chunks(gtf, threads=1000)
assert len(chunks) != 1000, 'get_chunks() failed to catch too many threads'

cpu = cpu_count()
chunks = get_chunks(gtf)
assert len(chunks) <= cpu, 'get_chunks() failed with no thread option provided for test GTF'

chunks = get_chunks(gtf, threads=2)
assert chunks[0].shape == (28,9) and chunks[1].shape == (36,9), 'get_chunks() failed output proper size chunks'


### Functional tests
from xpresspipe.gtfModify import longest_transcripts, protein_gtf, edit_gtf
gtf = pd.read_csv(
    str(__path__) + 'references/other/gtf_test_large.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)
gtf_edit_truncated = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=True,
    truncate_reference=True,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
truncate_truth = [
[1,	'ensembl',	'transcript',	69055,	70108],
[1,	'ensembl',	'exon',	69055,	70108],
[1,	'ensembl',	'CDS',	69136,	69990],
[1,	'ensembl',	'start_codon',	69091,	69093],
[1,	'ensembl',	'stop_codon',	70006,	70008],
[1,	'ensembl',	'five_prime_utr',	69055,	69090],
[1,	'ensembl',	'three_prime_utr',	70009,	70108],
[1,	'ensembl_havana',	'transcript',	450703,	451697],
[1,	'ensembl_havana',	'exon',	450703,	451697],
[1,	'ensembl_havana',	'CDS',	450758,	451633],
[1,	'ensembl_havana',	'start_codon',	451676,	451678],
[1,	'ensembl_havana',	'stop_codon',	450740,	450742],
[1,	'ensembl_havana',	'five_prime_utr',	451679,	451697],
[1,	'ensembl_havana',	'three_prime_utr',	450703,	450739]]
truncate_truth = pd.DataFrame(truncate_truth)
assert gtf_edit_truncated.iloc[:14,:5].equals(truncate_truth), 'edit_gtf() failed during truncation'
assert gtf_edit_truncated.shape == (10608, 9), 'Something went wrong during the functional test running through all GTF modifications'
