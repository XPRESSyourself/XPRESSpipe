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
[19,	'havana',	'transcript',	60951,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'exon',	        70928,	70976,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'exon',	        66346,	66499,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'exon',	        60951,	61894,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'transcript',	60951,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'exon',	        60951,	62000,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],
[19,	'havana',	'exon',	        69000,	70000,	'.',	'-',	'.',	'gene_id "ENSG00000282458";' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'+',	'.',	'gene_id "ENSG00000282495";' ],

[19,	'havana',	'gene',	        60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'transcript',	60951,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'exon',	        60951,	60960,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'exon',	        60970,	61500,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],
[19,	'havana',	'exon',	        62000,	63000,	'.',	'-',	'.',	'gene_id "ENSG00000282496";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'exon',	        60951,	61000,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],
[19,	'havana',	'exon',	        61005,	61010,	'.',	'-',	'.',	'gene_id "ENSG00000282498";' ],

[19,	'havana',	'gene',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'transcript',	60951,	61050,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'exon',	        60951,	61000,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
[19,	'havana',	'exon',	        61005,	61010,	'.',	'+',	'.',	'gene_id "ENSG00000282499";' ],
]
gtf = pd.DataFrame(gtf)

from xpresspipe.gtfTruncate import truncate_gtf
truncate_gtf(gtf)
