# Read Distribution

# Metagene
"""Test first step
"""
import os
import pandas as pd
from xpresspipe.complexity import make_complexity
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

# Complexity
args_dict = {}
args_dict['complexity'] = __path__ + 'complexity/'
args_dict['input'] = __path__ + 'complexity/input/'
args_dict['gtf'] = __path__ + 'references/other/gtf_test_large.gtf'
args_dict['output'] = __path__ + 'complexity/output/'
args_dict['log'] = ' '
args_dict['type'] = 'SE'
args_dict['experiment'] = 'test'
args_dict['path'] = str(__path__).replace('tests/', 'xpresspipe/')

make_complexity(args_dict)

assert os.path.isfile(args_dict['output'] + 'complexity/test_library_complexity_1_summary.pdf'), 'Complexity test failed to generate summary plot'
assert os.path.isfile(args_dict['output'] + 'complexity/metrics/sample.sort_dedupMarked_metrics.txt'), 'Complexity test failed to generate metrics file'

df = pd.read_csv(
    os.path.join(args_dict['output'], 'complexity/metrics/sample.sort_dedupMarked_metrics.txt'),
    sep='\t',
    index_col=0
)

""" 
                ID  geneLength  allCountsMulti  filteredCountsMulti  dupRateMulti  dupsPerIdMulti  RPKMulti  PKMMulti  allCounts  filteredCounts  dupRate  dupsPerId  RPK  RPKM
1  ENSG00000223972        1735               0                    0           NaN               0       0.0       0.0          0               0      NaN          0  0.0   0.0
2  ENSG00000227232        1351               0                    0           NaN               0       0.0       0.0          0               0      NaN          0  0.0   0.0
3  ENSG00000278267          68               0                    0           NaN               0       0.0       0.0          0               0      NaN          0  0.0   0.0
4  ENSG00000243485        1021               0                    0           NaN               0       0.0       0.0          0               0      NaN          0  0.0   0.0
5  ENSG00000284332         138               0                    0           NaN               0       0.0       0.0          0               0      NaN          0  0.0   0.0
"""

# Check ID column values 
id_values = ['ENSG00000223972', 'ENSG00000227232', 'ENSG00000278267', 'ENSG00000243485', 'ENSG00000284332']
df_genelens = [1735, 1351, 68, 1021, 138]

assert df.head()['ID'].tolist() == id_values, 'ID values are incorrect'
assert df.head()['geneLength'].tolist() == df_genelens, 'Gene lengths are incorrect'


"""
"60"	"ENSG00000187634"	4173	1	1	0	0	0.239635753654445	64.7839290766275	1	1	0	0	0.239635753654445	64.7839290766275
"61"	"ENSG00000188976"	5540	2	2	0	0	0.36101083032491	97.5968722154392	2	2	0	0	0.36101083032491	97.5968722154392
"""
# Check RPKM values for rows 60 and 61
rpk_multi_values = [0.239635753654445, 0.36101083032491]
rpkm_values = [64.7839290766275, 97.5968722154392]

assert df.iloc[59]['RPKMulti'] == rpk_multi_values[0], 'RPKMulti value is incorrect'
assert df.iloc[59]['RPKM'] == rpkm_values[0], 'RPKM value is incorrect'
assert df.iloc[60]['RPKMulti'] == rpk_multi_values[1], 'RPKMulti value is incorrect'
assert df.iloc[60]['RPKM'] == rpkm_values[1], 'RPKM value is incorrect'
