import os
import sys

__path__  =  os.path.dirname(os.path.realpath(__file__))
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'

# Test software install
assert os.system('xpresspipe --help') == 0, 'xpresspipe help call not working'
assert os.system('xpresspipe seRNAseq --help') == 0, 'seRNAseq help call not working'
assert os.system('xpresspipe peRNAseq --help') == 0, 'peRNAseq help call not working'
assert os.system('xpresspipe riboprof --help') == 0, 'riboprof help call not working'
assert os.system('xpresspipe trim --help') == 0, 'trim help call not working'
assert os.system('xpresspipe align --help') == 0, 'align help call not working'
assert os.system('xpresspipe count --help') == 0, 'count help call not working'
assert os.system('xpresspipe normalizeMatrix --help') == 0, 'normalizeMatrix help call not working'
assert os.system('xpresspipe diffxpress --help') == 0, 'diffxpress help call not working'
assert os.system('xpresspipe metagene --help') == 0, 'metagene help call not working'
assert os.system('xpresspipe readDistribution --help') == 0, 'readDistribution help call not working'
assert os.system('xpresspipe periodicity --help') == 0, 'periodicity help call not working'
assert os.system('xpresspipe complexity --help') == 0, 'complexity help call not working'
assert os.system('xpresspipe curateReference --help') == 0, 'curateReference help call not working'
assert os.system('xpresspipe makeReference --help') == 0, 'makeReference help call not working'
assert os.system('xpresspipe modifyGTF --help') == 0, 'modifyGTF help call not working'
assert os.system('xpresspipe rrnaProbe --help') == 0, 'rrnaProbe help call not working'
assert os.system('xpresspipe convertNames --help') == 0, 'convertNames help call not working'

# Test riboprof pipeline on test data
input = str(__path__) + 'riboprof_test'
output = str(__path__) + 'riboprof_out'
reference = str(__path__) + 'se_reference'
gtf = str(__path__) + 'se_reference/transcripts_longestTranscripts_proteinCoding_truncated.gtf'
os.system(
    'xpresspipe riboprof' \
    + ' -i ' + str(input) \
    + ' -o ' + str(output) \
    + ' -r ' + str(reference) \
    + ' -g ' + str(gtf) \
    + ' -e riboprof_test' \
    + ' -a CTGTAGGCACCATCAAT' \
    + ' --method RPKM' \
    + ' --sjdbOverhang 49')
