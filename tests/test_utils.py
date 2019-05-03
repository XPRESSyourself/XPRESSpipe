# Initialize functions
from xpresspipe.utils import check_directories, add_directory, get_files, get_probe_files, unzip_files, get_fasta
import os
import sys
__path__, xpresspipe_arguments  =  os.path.split(__file__)
__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

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

# Add directory
args_dict = {}
args_dict['log'] = ''
args_dict['old_dir'] = __path__
parent1 = 'old_dir'
parent2 = 'log'
child1 = 'test_make'
child2 = 'filename.txt'
child3 = 'test_make/'

assert os.path.exists(add_directory(args_dict, parent1, child1)[child1]) == True, 'add_directory() failed to add directory' # Make sure a properly formatted directory can be added
os.system('rm -r ' + str(args_dict[child1]))

try:
    add_directory(args_dict, parent2, child1) # Make sure Exception is caught when invalid parent directory is provided
except Exception:
    pass
else:
    raise Error('add_directory() failed')

try:
    add_directory(args_dict, parent1, child2) # Make sure Exception is caught when invalid child directory is provided
except Exception:
    pass
else:
    raise Error('add_directory() failed')

assert add_directory(args_dict, parent1, child3)['test_make'][-5:] == 'make/', 'add_directory() failed to properly format new directory name' # Make sure removes trailing '/' from new name
os.system('rm -r ' + str(args_dict[child1]))

# Get files




# Get prober files





# Unzip files



# Get FASTA files
files = []
def make_file(path, name):

    file = open(str(path) + str(name), 'w')
    file.close()
    files.append(str(path) + str(name))

make_file(__path__, 'fast1.txt')




for f in files:
    os.system('rm ' + str(f))
