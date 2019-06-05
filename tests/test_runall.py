import os
import sys

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'

os.system(
'python ' + str(__path__) + 'test_analysis.py'
)
os.system(
'python ' + str(__path__) + 'test_arguments.py'
)
os.system(
'python ' + str(__path__) + 'test_bam.py'
)
os.system(
'python ' + str(__path__) + 'test_flattenGTF.py'
)
os.system(
'python ' + str(__path__) + 'test_gtf.py'
)
os.system(
'python ' + str(__path__) + 'test_help.py'
)
os.system(
'python ' + str(__path__) + 'test_modifyGTF.py'
)
os.system(
'python ' + str(__path__) + 'test_modules.py'
)
os.system(
'python ' + str(__path__) + 'test_pipelines.py'
)
os.system(
'python ' + str(__path__) + 'test_truncateGTF.py'
)
os.system(
'python ' + str(__path__) + 'test_utils.py'
)
