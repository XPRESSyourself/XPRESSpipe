# Read Distribution

# Metagene
"""Test first step
"""
import os
import pandas as pd
from xpresspipe.metagene import run_metagene
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

args_dict = {
    'input': __path__ + 'data/',
    'output': __path__,
    'feature_type': 'exon',
    'experiment': 'test',
    'gtf': __path__ + 'references/other/gtf_metagene_test.gtf',
    'metagene': __path__ + 'metagene/',
    'metrics': __path__ + 'metagene/metrics/',
    'log': ' >> ' + __path__ + 'metagene/test.log 2>&1',
    'path': str(__path__).replace('tests/', 'xpresspipe/')
}
file = 'sample_small.bam'

os.system('mkdir ' + args_dict['metagene'])
os.system('mkdir ' + args_dict['metagene'] + 'metrics')

args = [file, args_dict]
run_metagene(args)

table = pd.read_csv(
    args_dict['metrics'] + 'sample_small.metaposit',
    sep='\t',
    index_col=0
)
table.sort_values(by=['meta_position'], inplace=True)
table = table.reset_index(drop=True)

assert table.seqnames.unique()[0] == 'ENST00000434826', 'Metagene test failed to find metaposition'
assert table.meta_position.tolist() == [
    103,
    308,
    493,
    602,
    933,
    944,
    1629,
    2012,
    2162,
    2427], 'Metagene test failed to find metaposition'

"""Test finisher
"""
from xpresspipe.metagene import finish_metagene
args = [file, args_dict]
finish_metagene(args, remove_outliers=False)

meta = pd.read_csv(
    args_dict['metrics'] + 'sample_small_metrics.txt',
    sep='\t'
)

assert meta.dropna()['representative transcript'].tolist() == [5, 13, 21, 26, 39, 40, 68, 84, 91]
assert meta.dropna()['metacount'].tolist() == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


args_dict['feature_type'] = 'CDS'
args = [file, args_dict]
run_metagene(args)
finish_metagene(args, remove_outliers=False)

meta = pd.read_csv(
    args_dict['metrics'] + 'sample_small_metrics.txt',
    sep='\t'
)

assert meta.dropna()['representative transcript'].tolist() == [6, 12, 30, 31, 69, 90, 98]
assert meta.dropna()['metacount'].tolist() == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

os.system('rm -r ' + args_dict['metagene'])

# Periodicity

# Complexity
