# Initialize functions
import os
import sys
import multiprocessing
import unittest
import xpresspipe as xp
__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

"""
check_inputs() -- General tests
"""
# Combination 0
test0_dict = {
    'max_processors': 1000000,
    'cmd': 'riboseq'} # Where input is more than available cores

try:
    t0 = xp.check_inputs(test0_dict)
except:
    raise Exception('check_inputs() failed checking max_processors')

assert type(t0['max_processors']) is int, 'check_inputs() failed checking max_processors' # Check cases where max_processors not None for int output

# Combination 1
test1_dict = {
    'max_processors': 'None',
    'cmd': 'riboseq',
    'output': __path__} # Where no limit is specified

try:
    t1 = xp.check_inputs(test1_dict)
except:
    raise Exception('check_inputs() failed checking max_processors')

# Combination 2
test2_dict = {
    'max_processors': multiprocessing.cpu_count(),
    'cmd': 'riboseq',
    'output': __path__} # Where limit equals available cores

try:
    t2 = xp.check_inputs(test2_dict)
except:
    raise Exception('check_inputs() failed checking max_processors')

assert type(t2['max_processors']) is int, 'check_inputs() failed checking max_processors'

# Combination 3
test3_dict = {
    'max_processors': multiprocessing.cpu_count() - 1,
    'cmd': 'riboseq',
    'output': __path__} # Where limit is less than available cores

try:
    t3 = xp.check_inputs(test3_dict)
except:
    raise Exception('check_inputs() failed checking max_processors')

assert type(t3['max_processors']) is int, 'check_inputs() failed checking max_processors'

"""
check_inputs() -- Directory formatting tests
"""
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

t = xp.check_inputs(test_dict)

assert truth_dict['input'] == t['input'], 'check_inputs() failed at directory formatting'
assert truth_dict['output'] == t['output'], 'check_inputs() failed at directory formatting'
assert truth_dict['reference'] == t['reference'], 'check_inputs() failed at directory formatting'

"""
check_inputs() -- adapter input tests
"""
# Test 0
adap1 = 'CTGTAGGCACCATCAAT'

test_dict = {
    'adapters': adap1,
    'cmd': 'riboseq'}

try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to catch non-list input')

# Test 1
adap2 = 'CTGTAGGCACCATCAAT CTGTAGGCACCATCAAG'
test_dict = {
    'adapters': adap2,
    'cmd': 'riboseq'}
try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to catch non-list input')

# Test 2
adap3 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCACCATCAAG']
test_dict = {
    'adapters': adap3,
    'cmd': 'riboseq'} # Check that expected case works
t = xp.check_inputs(test_dict)
assert t['adapters'] == adap3

# Test 3
adap4 = ['None']
test_dict = {
    'adapters': adap4,
    'cmd': 'riboseq'} # Check the None input works
t = xp.check_inputs(test_dict)
assert t['adapters'] == ['NONE']

# Test 4
adap5 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCACCATCAAG', 'ACCATCAAG']
test_dict = {
    'adapters': adap5,
    'cmd': 'riboseq'}
try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to catch input list with greater than two adapters')

# Test 5
adap6 = ['CTGTAGGCACCATCAAT', 'CTGTAGGCA262343CCATCAaAG']
test_dict = {
    'adapters': adap6,
    'cmd': 'riboseq'}
try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to make sure input adapters only contain valid characters')

# Test 6
adap7 = ['GCTCGCGCHATC']
test_dict = {
    'adapters': adap7,
    'cmd': 'riboseq'}
try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to make sure input adapters only contain valid characters')

# Test 7
adap8 = False
test_dict = {
    'adapters': adap8,
    'cmd': 'riboseq'}
try:
    t = xp.check_inputs(test_dict)
except Exception:
    pass
else:
    raise Exception('check_inputs() failed to make sure input adapters only contain valid characters')

"""
check_inputs() -- Log file formatting tests
"""
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

t0 = xp.check_inputs(test0_dict)
t1 = xp.check_inputs(test1_dict)

assert t0['output'] == t0['log_loc'], 'check_inputs() failed at log file formatting'
assert 'h_' in t0['log_file'], 'check_inputs() failed at log file formatting'
assert 'test' in t1['log_file'], 'check_inputs() failed at log file formatting'
