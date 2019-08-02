# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
import xpresspipe as xp
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

"""Test gtfFlatten functions"""

# Input file will test chunking to make sure GTF structure is conserved
gtf = str(__path__) + 'references/other/gtf_test_large.gtf'

"""
make_features() -- Test feature building
"""
# Test that dataframe is parsed for data correctly
gtf_frame = pd.read_csv(gtf,sep='\t',header=None, comment = '#',low_memory=False)

output = str(__path__) + 'references/other/gtf_test_large.fts'
gtf_feats = xp.make_features(gtf_frame, output)

# Test that data is output to expected position
try:
    # Test that output data was formatted correctly
    feats = pd.read_csv(output, sep='\t')
    assert feats.columns.tolist() == ['transcript','feature','length'], 'feature extraction failed'
    assert feats.shape == (42762, 3), 'feature extraction failed'
except:
    raise Exception('No output file detected for make_features()')

"""
make_index() -- Test index building
"""
from xpresspipe.buildIndex import make_index




"""
index_gtf() -- Test the main function
"""
from xpresspipe.buildIndex import index_gtf





















col = ['gene', 'strand', 'chromosome', 'start', 'end', 'coordinates', 'length']
gtf_truth = [['ENSG00000186092', '+', '1', 69055, 70108, [[69055, 70108]], 1054]]
gtf_truth = pd.DataFrame(gtf_truth, columns = col)
flat_gtf = flatten_reference(gtf)
assert flat_gtf.iloc[:1].equals(gtf_truth), (print(flat_gtf), print(gtf_truth), 'flatten_reference() failed with gtf dataframe input')

# Indexing truth set
flat = [['ENSG00000176695', '+', '19', 107104, 117102, [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]], 790],
        ['ENSG00000176696', '-', '19', 107104, 117102, [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]], 790],
        ['ENSG00000176697', '+', '19', 107104, 117102, [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]], 790],
        ['ENSG00000176698', '-', '19', 107104, 117102, [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]], 790],
        ['ENSG00000176699', '+', '20', 107104, 117102, [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]], 790]]
flat = pd.DataFrame(flat, columns = col)
coordinate_truth = [[[107104, 117102, '+', [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]],790],
    [107104, 117102, '-', [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]],790],
    [107104, 117102, '+', [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]],790],
    [107104, 117102, '-', [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]],790]],
    [[107104, 117102, '+', [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]],790]]]

# create_chromosome_index
from xpresspipe.gtfFlatten import create_chromosome_index
assert create_chromosome_index(flat) == {'19':0,'20':1}, 'create_chromosome_index() failed'

# create_coordinate_index
from xpresspipe.gtfFlatten import create_coordinate_index
assert create_coordinate_index(flat) == coordinate_truth, 'create_coordinate_index() failed'

# get_coding_length
from xpresspipe.gtfFlatten import get_coding_length
length_truth = [[107104, 107157], [107473, 107555], [110625, 110681], [116507, 117102]]
length = 54 + 83 + 57 + 596
assert get_coding_length(flat['coordinates'][0]) == length, 'get_coding_length() failed'


# Flat list
from xpresspipe.gtfFlatten import flat_list
L = ['asdf','asdfasdf','jionndwa','asnoaljca','asdkjha','pkcansasd']
L1 = [['asdf','asdfasdf'],['jionndwa','asnoaljca'],['asdkjha','pkcansasd']]
L2 = [['asdf','asdfasdf'],['jionndwa','asnoaljca'],['asdkjha']]
assert flat_list(L1) == L, 'flat_list() failed'
assert flat_list(L2) == L[:5], 'flat_list() failed'
