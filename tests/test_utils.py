# Initialize functions
from xpresspipe.utils import check_directories, add_directory, get_files, get_probe_files, unzip_files
import os
import sys
__path__, xpresspipe_arguments  =  os.path.split(__file__)
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests'
__path__ = __path__ + '/'

def make_file(path, name, file_list):

    file = open(str(path) + str(name), 'w')
    file.close()
    files.append(str(path) + str(name))

# Check directory
str1 = '/path/to/files'
str2 = '/path/to/files/'
str3 = './path'
str4 = './path/'
str5 = '/path/to/files/filename.txt'

assert check_directories(str1) == str2, 'check_directories() failed to properly format a directory'
assert check_directories(str2) == str2, 'check_directories() failed to properly format a directory'
assert check_directories(str3) == str4, 'check_directories() failed to properly format a directory'
assert check_directories(str4) == str4, 'check_directories() failed to properly format a directory'

try:
    output5 = check_directories(str5)
except:
    pass



# Get prober files



# Unzip files



# Get FASTA files
