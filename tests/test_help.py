import subprocess

from xpresspipe.test import run_command

# Test software install
run_command('xpresspipe --help')
run_command('xpresspipe seRNAseq --help')
run_command('xpresspipe peRNAseq --help')
run_command('xpresspipe riboseq --help')
run_command('xpresspipe trim --help')
run_command('xpresspipe align --help')
run_command('xpresspipe count --help')
run_command('xpresspipe normalizeMatrix --help')
run_command('xpresspipe diffxpress --help')
run_command('xpresspipe metagene --help')
run_command('xpresspipe readDistribution --help')
run_command('xpresspipe periodicity --help')
run_command('xpresspipe complexity --help')
run_command('xpresspipe curateReference --help')
run_command('xpresspipe makeReference --help')
run_command('xpresspipe modifyGTF --help')
run_command('xpresspipe rrnaProbe --help')
run_command('xpresspipe convertNames --help')
