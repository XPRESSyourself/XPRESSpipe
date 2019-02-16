"""
XPRESSpipe
An alignment and analysis pipeline for RNAseq data
alias: xpresspipe

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""



















"""
DESCRIPTION: Truncate user provided GTF transcript reference file

USAGE: Output file helpful for counting Ribosome Profiling reads to avoid calculating 5' bias

ASSUMPTIONS:
"""
def run_truncate(args_dict):


"""
DESCRIPTION: Compile counts tables from HTseq output files

ASSUMPTIONS:
FastQC has already been run on all files
Output folders from FastQC are in their own folder (path will be provided by user)
"""
def run_rrnaprobe(args_dict):
