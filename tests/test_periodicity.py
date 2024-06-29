# Read Distribution

# Metagene
"""Test first step
"""
import os
import pandas as pd
from xpresspipe.periodicity import make_periodicity
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

args_dict = {
    'input': __path__ + 'data/',
    'output': __path__,
    'feature_type': 'exon',
    'experiment': 'test',
    'bam_suffix': 'toTranscriptome.out.bam',
    'log': ' ',
    'path': str(__path__).replace('tests/', 'xpresspipe/')
}

# Periodicity
gtf = 'periodicity/references/filtered_transcripts.gtf.gz'
if os.path.isfile(__path__ + gtf):
    os.system('gunzip ' + __path__ + gtf)

args_dict['periodicity'] = __path__ + 'periodicity/'
args_dict['input'] = __path__ + 'periodicity/output/alignments_transcriptome/'
args_dict['gtf'] = __path__ + 'periodicity/references/filtered_transcripts.gtf'
args_dict['cdna_fasta'] = __path__ + 'periodicity/references/filtered_cds.fa'
args_dict['min_length'] = 17
args_dict['max_length'] = 0

make_periodicity(args_dict)

# Gzip the gtf file
os.system('gzip ' + args_dict['gtf'])

df = pd.read_csv(
    os.path.join(args_dict['output'], 'p_site_qc/metrics/SRR1795425_Aligned_small_metrics.txt'),
    sep='\t',
    index_col=0
).head()

"""
""	"transcript"	"psite"	"length"
"1"	"ENST00000001146"	673	30
"2"	"ENST00000001146"	3171	28
"3"	"ENST00000002125"	-242	31
"4"	"ENST00000002125"	-78	29
"5"	"ENST00000002125"	-62	32
"""
df_transcript_col = ['ENST00000001146', 'ENST00000001146', 'ENST00000002125', 'ENST00000002125', 'ENST00000002125']
df_psite_col = [673, 3171, -242, -78, -62]
df_length_col = [30, 28, 31, 29, 32]

assert df['transcript'].tolist() == df_transcript_col, 'Periodicity test failed to find correct transcripts'
assert df['psite'].tolist() == df_psite_col, 'Periodicity test failed to find correct psites'
assert df['length'].tolist() == df_length_col, 'Periodicity test failed to find correct lengths'

assert os.path.isfile(os.path.join(__path__, 'p_site_qc/codon_usage/test_codon_usage_1_summary.pdf')), 'Periodicity test failed to find codon usage summary'
assert os.path.isfile(os.path.join(__path__, 'p_site_qc/periodicity/test_periodicity_1_summary.pdf')), 'Periodicity test failed to find periodicity summary'