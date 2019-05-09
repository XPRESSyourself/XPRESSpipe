# Initialize functions
from xpresspipe.gtfModify import convert_gtf
import os
import sys
__path__, xpresspipe_arguments  =  os.path.split(__file__)
__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'


# Test conversion of GTF to fasta file for usage with Salmon
assert convert_gtf(str(__path__) + 'se_reference/transcripts.gtf', str(__path__) + 'se_reference/') == 0, 'convert_gtf() failed'
