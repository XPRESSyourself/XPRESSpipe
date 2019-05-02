# Initialize functions
from xpresspipe.arguments import check_inputs
import multiprocessing

# Test some input/output directory formatting checks
truth_dict = {
    'input': '/path/to/files/',
    'output': '/path/to/files2/',
    'reference': '/path/to/reference/'}

test_dict = {
    'input': '/path/to/files',
    'output': '/path/to/files2/',
    'reference': '/path/to/reference'}

t = check_inputs(test_dict)

assert truth_dict['input'] == t['input'], 'check_inputs() failed at directory formatting'
assert truth_dict['output'] == t['output'], 'check_inputs() failed at directory formatting'
assert truth_dict['reference'] == t['reference'], 'check_inputs() failed at directory formatting'

# Check max_processors argument
test0_dict = {'max_processors': 1000000} # Where input is more than available cores

test1_dict = {'max_processors': None} # Where no limit is specified

test2_dict = {'max_processors': multiprocessing.cpu_count()} # Where limit equals available cores

test3_dict = {'max_processors': multiprocessing.cpu_count() - 1} # Where limit is less than available cores

try:
    t0 = check_inputs(test0_dict)
    t1 = check_inputs(test1_dict)
    t2 = check_inputs(test2_dict)
    t3 = check_inputs(test3_dict)
except:
    raise Exception('check_inputs() failed checking max_processors')

assert type(t0['max_processors']) is int, 'check_inputs() failed checking max_processors' # Check cases where max_processors not None for int output
assert type(t2['max_processors']) is int, 'check_inputs() failed checking max_processors'
assert type(t3['max_processors']) is int, 'check_inputs() failed checking max_processors'

# Check adaptors
adap1 = 'CTGTAGGCACCATCAAT'
adap2 = 'CTGTAGGCACCATCAAT CTGTAGGCACCATCAAG'
adap3 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCACCATCAAG']
adap4 = None 


# Check log file formatting and location
test0_dict = {
    'input': '/path/to/files',
    'output': '/path/to/files2/',
    'reference': '/path/to/reference'}

test1_dict = {
    'input': '/path/to/files',
    'output': '/path/to/files2/',
    'reference': '/path/to/reference',
    'experiment': 'test'}

t0 = check_inputs(test0_dict)
t1 = check_inputs(test1_dict)

assert t0['output'] == t0['log_loc'], 'check_inputs() failed at log file formatting'
assert 'h_' in t0['log_file'], 'check_inputs() failed at log file formatting'
assert 'test' in t1['log_file'], 'check_inputs() failed at log file formatting'
