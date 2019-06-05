# Initialize functions
import os
import sys
import pandas as pd
__path__  =  os.path.dirname(os.path.realpath(__file__)) '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

# Read in BAM file to pandas dataframe
from xpresspipe.processBAM import read_bam
bam_file = str(__path__) + 'other/sample.bam'
bam_truth = [
['SRR1795425.22288830', 0, '1', 924081, 255, '32M', '*', 0, 0, 'CGGCGAGCCGGTCGTGGGACTGCCCCGGGCGC', 'CCFFFFFHHGHFHIGIJJGJJJJJJJJJHEDD', 'NH:i:1', 'HI:i:1', 'NM:i:1', 'MD:Z:11G20', 'AS:i:29'],
['SRR1795425.47844894', 16, '1', 952473, 255, '33M', '*', 0, 0, 'GCGTGCTGGTAGGCCACACCCGGCTCCAGGGCC', 'FAIGGEGGIHGIHFDC?CGEEHHHFHFDFFFC@', 'NH:i:1', 'HI:i:1', 'NM:i:0', 'MD:Z:33', 'AS:i:32'],
['SRR1795425.55229794', 16, '1', 957234, 255, '34M', '*', 0, 0, 'GAGCTGGTCTTTGTGCTCAGAGGCACGGCCTTTA', 'JJJJJJJJJJJJJJJJJJJJJJHHHHHFDFFFCC', 'NH:i:1', 'HI:i:1', 'NM:i:0', 'MD:Z:34', 'AS:i:33'],
['SRR1795425.47124117', 16, '1', 1217510, 255, '32M', '*', 0, 0, 'GCTCAAAACTCCTCGTGCACGCTGCGCGCGTA', 'IICGJJIHIHJJJJJJJJJJHHGHHFFFFFCC', 'NH:i:1', 'HI:i:1', 'NM:i:0', 'MD:Z:32', 'AS:i:31'],
['SRR1795425.63318546', 16, '1', 1312326, 255, '22M', '*', 0, 0, 'CCAGCTCTTTGAGGGCTTGCTC', 'JJIHHIHIHEHHHHHFFFFFCC', 'NH:i:1', 'HI:i:1', 'NM:i:0', 'MD:Z:22', 'AS:i:21'],
]
bam_truth = pd.DataFrame(bam_truth)

bam = read_bam(bam_file)
assert bam.head().equals(bam_truth), 'read_bam() failed'

bam = read_bam(bam_file, threads=2)
assert bam.head().equals(bam_truth), 'read_bam() failed with multiprocessing included'

# Sample bam file
from xpresspipe.processBAM import bam_sample
bam_small = bam_sample(bam, 3)
assert bam_small.shape == (3,16), 'bam_sample() failed '

bam_small = bam_sample(bam, 135)
assert bam_small.shape == (135,16), 'bam_sample() failed '
