import os
import sys
import resource
import gc

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

from xpresspipe.processBAM import read_bam
from xpresspipe.parallel import run_pools, run_pools_old
from xpresspipe.utils import get_files

def woop(bam):

    print('3Mem:'  + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    bam = bam[(bam[9].str.len() == 28)]
    return bam[9]

def woo(args):
    file, args_dict = args[0], args[1] # Parse args

    print('1starting: ' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    # Read in indexed bam file
    bam = read_bam(
        str(args_dict['input']) + str(file),
        threads=5)
    print('2file read' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

    bam_coordinates = woop(bam)

    print('4Cleaning:' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    bam = None # Some clean-up
    gc.collect()
    print('5Post:' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

args_dict = {}
args_dict['workers'] = 2
args_dict['input'] = '/Users/jordan/Desktop/xpressyourself_testing/7155666_export/alignments/aligned_merged/'

file_list = get_files(
    args_dict['input'],
    ['_merged.bam'])
args_iter = [[file, args_dict] for file in file_list]
print(file_list)
print('NEW---------------------')
run_pools(
    woo,
    args_iter,
    args_dict)

print('OLD---------------------')
run_pools_old(
    woo,
    args_iter,
    args_dict)
