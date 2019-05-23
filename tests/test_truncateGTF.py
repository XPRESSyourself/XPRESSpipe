from __future__ import print_function

import os
import sys
import pandas as pd
import numpy as np
__path__  =  os.path.dirname(os.path.realpath(__file__))
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'

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
gtf = pd.DataFrame(gtf)

gtf_truth0 = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'exon',	        70943,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'exon',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'exon',	        60951,	61849,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60966,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	69955,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'exon',	        61006,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62985,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60976,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
]
gtf_truth0 = pd.DataFrame(gtf_truth0)

gtf_truth1 = [
[19,	'havana',	'gene',	        60951,	71626,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458";' ],
[19,	'havana',	'exon',	        70933,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'exon',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'exon',	        60951,	61844,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60956,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        69000,	69950,	'.',	'-',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460";' ],
[19,	'havana',	'exon',	        61011,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST00000282460"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461";' ],
[19,	'havana',	'exon',	        60956,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62950,	'.',	'-',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST00000282461"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ]
]
gtf_truth1 = pd.DataFrame(gtf_truth1)

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
[19,	'havana',	'exon',	        70933,	70940,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "1"' ],
[19,	'havana',	'exon',	        70950,	70955,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "2"' ],
[19,	'havana',	'exon',	        70955,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282458"; exon_number "3"' ],

[19,	'havana',	'transcript',	60956,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60956,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "4"' ],
[19,	'havana',	'exon',	        60965,	60970,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "5"' ],
[19,	'havana',	'start',	    60956,	60959,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        60975,	60980,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "6"' ],
[19,	'havana',	'exon',	        60990,	64000,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "7"' ],
[19,	'havana',	'exon',	        64050,	64055,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "8"' ],
[19,	'havana',	'exon',	        64060,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "9"' ],

[19,	'havana',	'gene',	        60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260";' ],
[19,	'havana',	'exon',	        60000,	62000,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "2"' ],
[19,	'havana',	'exon',	        62970,	62971,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "3"' ],
[19,	'havana',	'exon',	        62975,	62980,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "4"' ],
[19,	'havana',	'exon',	        62990,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "5"' ],

[19,	'havana',	'gene',	        60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261";' ],
[19,	'havana',	'exon',	        60900,	60945,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "1"' ],
[19,	'havana',	'exon',	        60950,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262";' ],
[19,	'havana',	'exon',	        60900,	60915,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "1"' ],
[19,	'havana',	'exon',	        60930,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "2"' ],
[19,	'havana',	'exon',	        62000,	62045,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "3"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'transcript',	60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST0000028263";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282498"; transcript_id "ENST0000028263"; exon_number "1"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'transcript',	60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST0000028264";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282499"; transcript_id "ENST0000028264"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265";' ],
[19,	'havana',	'exon',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266";' ],
[19,	'havana',	'exon',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266"; exon_number "1"' ],

[19,	'havana',	'gene',	        60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505";' ],
[19,	'havana',	'transcript',	60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267";' ],
[19,	'havana',	'exon',	        60956,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "1"' ],
[19,	'havana',	'exon',	        60970,	62940,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "2"' ],
[19,	'havana',	'exon',	        62945,	62950,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "3"' ],
[19,	'havana',	'exon',	        62950,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "4"' ]
]
gtf2 = pd.DataFrame(gtf2)

gtf_truth2 = [
[19,	'havana',	'gene',	        60956,	70960,	'.',	'+',	'.',	'gene_id "ENSG00000282458";' ],

[19,	'havana',	'transcript',	60956,	64065,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'start',	    60956,	60959,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459";' ],
[19,	'havana',	'exon',	        61021,	63995,	'.',	'+',	'.',	'gene_id "ENSG00000282458"; transcript_id "ENST00000282459"; exon_number "7"' ],

[19,	'havana',	'gene',	        60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60000,	62995,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260";' ],
[19,	'havana',	'exon',	        60045,	61996,	'.',	'+',	'.',	'gene_id "ENSG00000282495"; transcript_id "ENST0000028260"; exon_number "2"' ],

[19,	'havana',	'gene',	        60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60900,	62015,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261";' ],
[19,	'havana',	'exon',	        60950,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282496"; transcript_id "ENST0000028261"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262";' ],
[19,	'havana',	'exon',	        60930,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282497"; transcript_id "ENST0000028262"; exon_number "2"' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265";' ],
[19,	'havana',	'exon',	        60915,	60916,	'.',	'-',	'.',	'gene_id "ENSG00000282500"; transcript_id "ENST0000028265"; exon_number "1"' ],

[19,	'havana',	'gene',	        60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501";' ],
[19,	'havana',	'transcript',	60900,	60961,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266";' ],
[19,	'havana',	'exon',	        60945,	60946,	'.',	'+',	'.',	'gene_id "ENSG00000282501"; transcript_id "ENST0000028266"; exon_number "1"' ],

[19,	'havana',	'gene',	        60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505";' ],
[19,	'havana',	'transcript',	60956,	62955,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267";' ],
[19,	'havana',	'exon',	        60981,	62905,	'.',	'-',	'.',	'gene_id "ENSG00000282505"; transcript_id "ENST0000028267"; exon_number "2"' ]
]
gtf_truth2 = pd.DataFrame(gtf_truth2)


from xpresspipe.gtfTruncate import truncate_gtf
# Test0
truncated_gtf0 = truncate_gtf(gtf)
truncated_gtf0 = truncated_gtf0.reset_index(drop=True)

#print(gtf)
#print('------------------')
#print(gtf_truth0)
#print('------------------')
#print(truncated_gtf0)

#print('???????????????/')

assert truncated_gtf0.equals(gtf_truth0), run_debugger(range(0,18), gtf_truth0, truncated_gtf0)

# Test1
truncated_gtf1 = truncate_gtf(gtf, _5prime=50, _3prime=5)
truncated_gtf1 = truncated_gtf1.reset_index(drop=True)

#print(gtf)
#print('----------Truth--------')
#print(gtf_truth1)
#print('---------Test---------')
#print(truncated_gtf1)

#print('???????????????/')

assert truncated_gtf1.equals(gtf_truth1), run_debugger(range(0,19), gtf_truth1, truncated_gtf1)

# Test2
truncated_gtf2 = truncate_gtf(gtf2)

#print(gtf2)
#print('----------Truth--------')
#print(gtf_truth2)
#print('---------Test---------')
#print(truncated_gtf2)

#print('???????????????/')

assert gtf_truth2.equals(truncated_gtf2), run_debugger(range(0,24), truncated_gtf2, gtf_truth2)
