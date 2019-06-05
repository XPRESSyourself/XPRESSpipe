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
    str(__path__) + 'other/gtf_test.gtf',
    sep = '\t',
    header = None,
    comment = '#',
    low_memory = False)

# Test gtfFlatten functions
from xpresspipe.gtfFlatten import flatten_reference

#df = flatten_reference(str(__path__) + 'se_reference/transcripts.gtf')
#df.head()


from xpresspipe.gtfFlatten import create_chromosome_index, create_coordinate_index, flat_list, get_coding_length, make_flatten, flatten_reference
