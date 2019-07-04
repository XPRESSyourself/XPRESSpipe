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
    str(__path__) + 'other/gtf_test.gtf',
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
assert gtf_edit.shape == (10608, 9), 'Something went wrong during the functional test running through all GTF modifications'

### Unit tests
# Get longest transcript records only for each gene
# Check CCDS proteins were prioritized over others
assert gtf_long.loc[(gtf_long[2] == 'transcript') & (gtf_long[8].str.contains('CCDS'))].shape == (398, 9), 'failed to catch CCDS domains'

# Do general spot check
gtf_longest_truth2 = [
[1,	'havana',	'transcript',	11869,	14409],
[1,	'havana',	'transcript',	14404,	29570],
[1,	'mirbase',	'transcript',	17369,	17436],
[1,	'havana',	'transcript',	29554,	31097],
[1,	'mirbase',	'transcript',	30366,	30503],
[1,	'havana',	'transcript',	34554,	36081]]
gtf_longest_truth2 = pd.DataFrame(gtf_longest_truth2, index=[0,4,16,18,22,24])
assert gtf_long.iloc[0:25,].loc[gtf_long[2] == 'transcript'].iloc[:,0:5].equals(gtf_longest_truth2), 'failed to pass general check of transcript choosing'

# Get protein coding only records
gtf_protein = protein_gtf(gtf)
assert gtf_protein.shape == (46427, 9), 'protein_gtf() failed'

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
assert gtf_edit.shape == (10694, 9), 'edit_gtf() failed'

gtf_edit = edit_gtf(
    gtf,
    longest_transcript=False,
    protein_coding=True,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
assert gtf_edit.shape == (46427,9), 'edit_gtf() failed during protein_coding only output'


gtf_edit = edit_gtf(
    gtf,
    longest_transcript=True,
    protein_coding=False,
    truncate_reference=False,
    _5prime=45,
    _3prime=15,
    output=False,
    threads=None)
gtf_longest_truth1 = [
[1,	'ensembl', 'transcript',	69055,	70108],
[1,	'ensembl_havana',	'transcript',	450703,	451697]]
gtf_longest_truth1 = pd.DataFrame(gtf_longest_truth1, index=[34,101])
assert gtf_edit.iloc[0:125,].loc[(gtf_edit[2] == 'transcript') & (gtf_edit[8].str.contains('CCDS'))].iloc[:,0:5].equals(gtf_longest_truth1), 'failed to catch CCDS domains'
assert gtf_edit.iloc[0:25,].loc[gtf_edit[2] == 'transcript'].iloc[:,0:5].equals(gtf_longest_truth2), 'failed to pass general check of transcript choosing'


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
