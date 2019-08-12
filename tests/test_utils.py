# Initialize functions
from xpresspipe.utils import check_directories, add_directory, get_files, unzip_files, get_fasta
import os
import sys
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'


def make_file(path, name, file_list):

    file = open(str(__path__) + str(name), 'w')
    file.close()
    files.append(str(__path__) + str(name))

# Check directory
str1 = __path__
str2 = __path__
str5 = __path__ + 'file.txt'

assert check_directories(str1) == str2, 'check_directories() failed to properly format a directory'

try:
    output5 = check_directories(str5)
except:
    pass



# Get prober files



# Unzip files



# Get FASTA files
fasta_path = str(__path__) + 'references/fasta_test/'
fasta_files = get_fasta(fasta_path)
fasta_truth = str(fasta_path) + 'genome_fasta/fast1.fasta ' + str(fasta_path) + 'genome_fasta/fast2.fasta ' + str(fasta_path) + 'genome_fasta/fast3.fasta'
assert fasta_files == fasta_truth, 'get_fasta() failed to recurse through children directories to find fasta files'

fasta_path2 = str(__path__) + 'references/fasta_test2/'
fasta_files2 = get_fasta(fasta_path2)
fasta_truth2 = str(fasta_path2) + 'fast1.fasta ' + str(fasta_path2) + 'fast2.fasta ' + str(fasta_path2) + 'fast3.fasta'
assert fasta_files2 == fasta_truth2, 'get_fasta() failed to recurse through children directories to find fasta files'
