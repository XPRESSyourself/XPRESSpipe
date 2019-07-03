# Initialize functions
import os
import sys
import pandas as pd
import numpy as np
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
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


from xpresspipe.geneCoverage import cov_coordinates
coordinates_truth = [
['1', 924081, 924113],
['1', 952473, 952506],
['1', 957234, 957268],
['1', 1217510, 1217542],
['1', 1312326, 1312348],
]

bam_coordinates = cov_coordinates(bam_truth)
assert np.array_equal(bam_coordinates, coordinates_truth), 'Failed to parse out relevant coverage bam stats'


from xpresspipe.geneCoverage import get_cov_position
truth_list = [210,40,4,5]

coordinates = [[35,41],[210,212],[0,10],[40120,40122]]
positions = [40,210,1320,40123,23,-1,4,5]
test_list = []
for position in positions:
    p = get_cov_position(position, coordinates)
    if p != None:
        test_list.append(p)

assert set(truth_list) == set(test_list)


from xpresspipe.geneCoverage import get_cov_profile
coordinate_index = [[[4490468, 4587469, '+', [[4490468, 4490770], [4544567, 4544707], [4561449, 4561541], [4564344, 4564458], [4566047, 4566089], [4567669, 4567767], [4572204, 4572388], [4573907, 4574014], [4576001, 4576123], [4576569, 4576763], [4583038, 4583172], [4585312, 4587469]], 3686]]]
aligned_reads_index = [
    ['1',4490470,4490490],
    ['1',4490450,4490466],
    ['1',4490470,4490500],
    ['1',4544707,4544709],
    ['1',4544708,4544710]]

metagene_profile = get_cov_profile(
    aligned_reads_index,
    coordinate_index)


meta_index = metagene_profile.index.tolist()
assert min(meta_index) == coordinate_index[0][0][0]
assert max(meta_index) == coordinate_index[0][0][1]
assert metagene_profile['raw_count'].sum() == 53
