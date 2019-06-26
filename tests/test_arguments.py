# Initialize functions
import os
import sys
import multiprocessing
import unittest
from xpresspipe.arguments import check_inputs
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

# Check max_processors argument
test0_dict = {
    'max_processors': 1000000,
    'cmd': 'riboseq'} # Where input is more than available cores
test1_dict = {
    'max_processors': 'None',
    'cmd': 'riboseq'} # Where no limit is specified
test2_dict = {
    'max_processors': multiprocessing.cpu_count(),
    'cmd': 'riboseq'} # Where limit equals available cores
test3_dict = {
    'max_processors': multiprocessing.cpu_count() - 1,
    'cmd': 'riboseq'} # Where limit is less than available cores

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

# Test some input/output directory formatting checks
truth_dict = {
    'input': __path__,
    'output': __path__,
    'reference': __path__,
    'cmd': 'riboseq'}
test_dict = {
    'input': __path__,
    'output': __path__,
    'reference': __path__,
    'cmd': 'riboseq'}

t = check_inputs(test_dict)

assert truth_dict['input'] == t['input'], 'check_inputs() failed at directory formatting'
assert truth_dict['output'] == t['output'], 'check_inputs() failed at directory formatting'
assert truth_dict['reference'] == t['reference'], 'check_inputs() failed at directory formatting'

# Check adaptors
adap1 = 'CTGTAGGCACCATCAAT'
adap2 = 'CTGTAGGCACCATCAAT CTGTAGGCACCATCAAG'
adap3 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCACCATCAAG']
adap4 = ['None']
adap5 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCACCATCAAG', 'ACCATCAAG']
adap6 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCA262343CCATCAaAG']
adap7 = ['GCTCGCGCHATC']
adap8 = False

test_dict = {
    'adaptors': adap1,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to catch non-list input')


test_dict = {
    'adaptors': adap2,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to catch non-list input')

test_dict = {
    'adaptors': adap3,
    'cmd': 'riboseq'} # Check that expected case works
t = check_inputs(test_dict)
assert t['adaptors'] == adap3

test_dict = {
    'adaptors': adap4,
    'cmd': 'riboseq'} # Check the None input works
t = check_inputs(test_dict)
assert t['adaptors'] == ['NONE']

test_dict = {
    'adaptors': adap5,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to catch input list with greater than two adaptors')

test_dict = {
    'adaptors': adap6,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to make sure input adaptors only contain valid characters')

test_dict = {
    'adaptors': adap7,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to make sure input adaptors only contain valid characters')

test_dict = {
    'adaptors': adap8,
    'cmd': 'riboseq'}
try:
    t = check_inputs(test_dict)
except Exception:
    pass
else:
    raise Error('check_inputs() failed to make sure input adaptors only contain valid characters')

# Check log file formatting and location
test0_dict = {
    'input': __path__,
    'output': __path__,
    'reference': __path__,
    'cmd': 'riboseq'}

test1_dict = {
    'input': __path__,
    'output': __path__,
    'reference': __path__,
    'experiment': 'test',
    'cmd': 'riboseq'}

t0 = check_inputs(test0_dict)
t1 = check_inputs(test1_dict)

assert t0['output'] == t0['log_loc'], 'check_inputs() failed at log file formatting'
assert 'h_' in t0['log_file'], 'check_inputs() failed at log file formatting'
assert 'test' in t1['log_file'], 'check_inputs() failed at log file formatting'
