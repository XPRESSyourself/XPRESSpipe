# Initialize functions
from xpresspipe.utils import check_directories, add_directory, get_files, unzip_files
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
