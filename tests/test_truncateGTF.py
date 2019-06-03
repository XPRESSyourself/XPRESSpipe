from __future__ import print_function

import os
import sys
import pandas as pd
import numpy as np
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

"""
# Test the following
- different lengths for ends
- transcript too short and is removed
- exon too short on either side and you have to continue to next exon

"""

def run_debugger(
    range,
    gtf1,
    gtf2):

    print('Error: truncate_gtf() failed')
    for x in range:
        if x < 10:
            print('0' + str(x),':',gtf1.iloc[x].equals(gtf2.iloc[x]))
        else:
            print(str(x),':',gtf1.iloc[x].equals(gtf2.iloc[x]))


"""truncateGTF functions"""
from xpresspipe.gtfTruncate import truncate_gtf
gtf_base = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'exon',	        70928,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'exon',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'exon',	        60951,	61894,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462";' ],
[19,	'havana',	'exon',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462"; exon_number "1"' ],
[19,	'havana',	'exon',	        61005,	61010,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463";' ],
[19,	'havana',	'exon',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463"; exon_number "1"' ],
[19,	'havana',	'exon',	        61005,	61010,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463"; exon_number "2"' ]
]
gtf_base = pd.DataFrame(gtf_base)

# Test0
truncated_gtf_base = truncate_gtf(gtf_base)
truncated_gtf_base = truncated_gtf_base.reset_index(drop=True)

assert truncated_gtf_base.equals(gtf_base), 'failed to skip record where none are CDS'


# Input file will test chunking to make sure GTF structure is conserved
# for 45 start and 15 stop....
# ENSG00000282458: Normal
# ENSG00000282495: start is too short for plus
# ENSG00000282496: end is too short only  for minus
# ENSG00000282498: Full record deleted because too short
# ENSG00000282499: Full record deleted because too short
gtf = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'CDS',	        70928,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'CDS',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'CDS',	        60951,	61894,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'CDS',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],
[19,	'havana',	'CDS',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'CDS',	        60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60970,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62000,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],
[19,	'havana',	'CDS',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462";' ],
[19,	'havana',	'CDS',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462"; exon_number "1"' ],
[19,	'havana',	'CDS',	        61005,	61010,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST00000282462"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463";' ],
[19,	'havana',	'CDS',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463"; exon_number "1"' ],
[19,	'havana',	'CDS',	        61005,	61010,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST00000282463"; exon_number "2"' ]
]
gtf = pd.DataFrame(gtf)

gtf_truth0 = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'CDS',	        70943,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'CDS',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'CDS',	        60951,	61849,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'CDS',	        60966,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],
[19,	'havana',	'CDS',	        69000,	69955,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'CDS',	        61006,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62000,	62985,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],
[19,	'havana',	'CDS',	        62006,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
]
gtf_truth0 = pd.DataFrame(gtf_truth0)

truncated_gtf0 = truncate_gtf(gtf, _5prime=45, _3prime=15)
truncated_gtf0 = truncated_gtf0.reset_index(drop=True)

assert truncated_gtf0.equals(gtf_truth0), run_debugger(range(0,19), gtf_truth0, truncated_gtf0)


gtf_truth1 = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'CDS',	        70933,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'CDS',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'CDS',	        60951,	61844,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'CDS',	        60956,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],
[19,	'havana',	'CDS',	        69000,	69950,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'CDS',	        61011,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60956,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],
[19,	'havana',	'CDS',	        62000,	62950,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ]
]
gtf_truth1 = pd.DataFrame(gtf_truth1)

truncated_gtf1 = truncate_gtf(gtf, _5prime=50, _3prime=5)
truncated_gtf1 = truncated_gtf1.reset_index(drop=True)

assert truncated_gtf1.equals(gtf_truth1), run_debugger(range(0,19), gtf_truth1, truncated_gtf1)


# Add more test cases where all exons get lost
#1. All get trimmed out
#2. 3 off front, two off back
#3. 3 off the back
#4. Front and back a perfect trim (+)
#5. Front and back a perfect trim (-)
#6. One exon and too short (-)
#7. One exon and too short (+)
#8. One exon and 1bp left (-)
#9. One exon and 1bp left (+)
#10. Trim one and remainder off top, 2 and remainder of bottom, leave some, and a (-) stranded transcript as last record
gtf2 = [
[19,	'havana',	'gene',	        60956,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	70933,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'CDS',	        70933,	70940,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'CDS',	        70950,	70955,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'CDS',	        70955,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60956,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'CDS',	        60956,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'CDS',	        60965,	60970,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],
[19,	'havana',	'start',	    60956,	60959,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'CDS',	        60975,	60980,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "6"' ],
[19,	'havana',	'CDS',	        60990,	64000,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "7"' ],
[19,	'havana',	'CDS',	        64050,	64055,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "8"' ],
[19,	'havana',	'CDS',	        64060,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "9"' ],

[19,	'havana',	'gene',	        60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260";' ],
[19,	'havana',	'CDS',	        60000,	62000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62970,	62971,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "3"' ],
[19,	'havana',	'CDS',	        62975,	62980,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "4"' ],
[19,	'havana',	'CDS',	        62990,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "5"' ],

[19,	'havana',	'gene',	        60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261";' ],
[19,	'havana',	'CDS',	        60900,	60945,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60950,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62000,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262";' ],
[19,	'havana',	'CDS',	        60900,	60915,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60930,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62000,	62045,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'transcript',	60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST0000028263";' ],
[19,	'havana',	'CDS',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST0000028263"; exon_number "1"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'transcript',	60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST0000028264";' ],
[19,	'havana',	'CDS',	        60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST0000028264"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265";' ],
[19,	'havana',	'CDS',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266";' ],
[19,	'havana',	'CDS',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266"; exon_number "1"' ],

[19,	'havana',	'gene',	        60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505";' ],
[19,	'havana',	'transcript',	60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267";' ],
[19,	'havana',	'CDS',	        60956,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "1"' ],
[19,	'havana',	'CDS',	        60970,	62940,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "2"' ],
[19,	'havana',	'CDS',	        62945,	62950,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "3"' ],
[19,	'havana',	'CDS',	        62950,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "4"' ]
]
gtf2 = pd.DataFrame(gtf2)

gtf_truth2 = [
[19,	'havana',	'gene',	        60956,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458";' ],

[19,	'havana',	'transcript',	60956,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'start',	    60956,	60959,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'CDS',	        61021,	63995,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "7"' ],

[19,	'havana',	'gene',	        60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260";' ],
[19,	'havana',	'CDS',	        60045,	61996,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "2"' ],

[19,	'havana',	'gene',	        60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261";' ],
[19,	'havana',	'CDS',	        60950,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262";' ],
[19,	'havana',	'CDS',	        60930,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265";' ],
[19,	'havana',	'CDS',	        60915,	60916,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266";' ],
[19,	'havana',	'CDS',	        60945,	60946,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266"; exon_number "1"' ],

[19,	'havana',	'gene',	        60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505";' ],
[19,	'havana',	'transcript',	60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267";' ],
[19,	'havana',	'CDS',	        60981,	62905,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "2"' ]
]
gtf_truth2 = pd.DataFrame(gtf_truth2)

truncated_gtf2 = truncate_gtf(gtf2)
truncated_gtf2 = truncated_gtf2.reset_index(drop=True)

assert gtf_truth2.equals(truncated_gtf2), run_debugger(range(0,24), truncated_gtf2, gtf_truth2)
