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
IMPORT DEPENDENCIES
"""
import os, sys
import argparse
import multiprocessing
from textwrap import dedent
from .utils import check_directories

"""
INITIALIZATION PARAMETERS
"""
#Retrieve path for scripts used in this pipeline, appended to argument dictionary for every function
__path__, xpresspipe_arguments = os.path.split(__file__)

#Set default values for arguments
DEFAULT_READ_MIN = 18
DEFAULT_READ_QUALITY = 28
DEFAULT_MAX_PROCESSORS = None

descrip = """\
The XPRESSpipe submodules can be accessed by executing:
    'xpresspipe module_name'


"""

"""
DESCRIPTION: Check if file is comma separated
"""
def check_csv(file, sep=','):

    #Make sure file is a csv (for use with converting count tables to gene names and RPKM values)
    if str(file).endswith('.csv'):
        pass
    else:
        raise argparse.ArgumentTypeError('Input file type incorrect')
        sys.exit(1)

    return str(file)

"""
DESCRIPTION: Check arguments provided by user
"""
def check_inputs(args_dict):

    #Check user-provided directory formatting
    if 'input' in args_dict:
        args_dict['input'] = check_directories(args_dict['input'])
    if 'output' in args_dict:
        args_dict['output'] = check_directories(args_dict['output'])
    if 'reference in args_dict':
        args_dict['reference'] = check_directories(args_dict['reference'])

    #Check max_processor input
    if 'max_processors' in args_dict and args_dict['max_processors'] != None:
        args_dict['max_processors'] = int(args_dict['max_processors'])

        if multiprocessing.cpu_count() < args_dict['max_processors']:
            raise Exception('Cannot specify more cores than are available -- Specified ' + str(args_dict['max_processors']) + ' cores, only ' + str(multiprocessing.cpu_count()) + ' available')

    #Check number of adaptors provided
    if 'adaptors' in args_dict:
        if type(args_dict['adaptors']) != list:
            raise Exception('Adaptors must be provided as a list of strings')

        for x in args_dict['adaptors']:
            if type(x) != str:
                raise Exception('Adaptors must be provided as a list of strings')

        if len(args_dict['adaptors'] > 2:
            raise Exception('A maximum of 2 adaptors may be provided')

"""
DESCRIPTION: Get user arguments to determine sub-module to run and arguments provided
"""
def get_arguments(args, __version__):

    if args is None:
        args = sys.argv[1:] #requires user input

    """
    INITIALIZE PARSER
    """
    parser = argparse.ArgumentParser(prog='RiboPipe', description=dedent(descrip), formatter_class=argparse.RawDescriptionHelpFormatter)
    #optional arguments
    parser.add_argument(
        '-v', '--version',
        help='Print installed version to stout',
        action='version',
        version='%(prog)s ' + str(__version__)
        )

    """
    MODULE SUBPARSER PROGRAMS
    """
    subparser = parser.add_subparsers(title='Sub-modules', description='Choose one of the following:', dest='cmd')

    #TRIM subparser
    trim_parser = subparser.add_parser('trim', description='Trim RNAseq reads of adaptors and for quality', add_help=False)
    #Required arguments
    trim_reqs = trim_parser.add_argument_group('required arguments')
    trim_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory -- if paired-end, file names should be exactly the same except for r1/r2.fastq or similar suffix',
        required=True)
    trim_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    #Optional arguments
    trim_opts = trim_parser.add_argument_group('optional arguments')
    trim_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    trim_opts.add_argument(
        '-a', '--adaptors',
        metavar='<list>',
        nargs='+',
        default=None,
        help='Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if none are provided, software will attempt to auto-detect adaptors -- if \"POLYX\" is provided as a single string in the list, polyX adaptors will be trimmed',
        required=False)
    trim_opts.add_argument(
        '-q', '--quality',
        metavar='<int>',
        default=DEFAULT_READ_QUALITY,
        help='PHRED read quality threshold (default: %s)' % DEFAULT_READ_QUALITY,
        required=False)
    trim_opts.add_argument(
        '--min_length',
        metavar='<int>',
        default=DEFAULT_READ_MIN,
        help='Minimum read length threshold to keep for reads (default: %s)' % DEFAULT_READ_MIN,
        required=False)
    trim_opts.add_argument(
        '-m', '--max_processors',
        help='Number of max processors to use for tasks (default: No limit)',
        metavar='<int>',
        default=DEFAULT_MAX_PROCESSORS,
        required=False)

    #ALIGN subparser
    align_parser = subparser.add_parser('align', description='Align RNAseq reads to reference genome', add_help=False)
    #Required arguments
    align_reqs = align_parser.add_argument_group('required arguments')
    align_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory -- if paired-end, file names should be exactly the same except for r1/r2.fastq or similar suffix',
        required=True)
    align_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    align_reqs.add_argument(
        '-t', '--type',
        help='Sequencing type (\"SE\" for single-end, \"PE\" for paired-end)',
        required=True)
    align_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        required=True)
    #Optional arguments
    align_opts = align_parser.add_argument_group('optional arguments')
    align_opts.add_argument(
        '--output_bed',
        help='Include flag to output BED files for each aligned file',
        action='store_true',
        required=False)
    align_opts.add_argument(
        '--output_bigwig',
        help='Include flag to output bigwig files for each aligned file',
        action='store_true',
        required=False)

    #COUNT subparser
    _parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    _reqs = _parser.add_argument_group('required arguments')
    _opts.add_argument()
    #Optional arguments
    _opts = _parser.add_argument_group('optional arguments')
    _opts.add_argument(
        '--truncate',
        help='Count reads to truncated reference (truncated reference should be added to the \"reference\" parent directory)',
        type=bool,
        required=False)
    _opts.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=False)

    #QUALITY subparser
    _parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    _reqs = _parser.add_argument_group('required arguments')
    _opts.add_argument()
    #Optional arguments
    _opts = _parser.add_argument_group('optional arguments')
    _opts.add_argument()

    #TRUNCATE subparser
    _parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    _reqs = _parser.add_argument_group('required arguments')
    _opts.add_argument()
    #Optional arguments
    _opts = _parser.add_argument_group('optional arguments')
    _opts.add_argument()

    #RRNAPROBE subparser
    _parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    _reqs = _parser.add_argument_group('required arguments')
    _opts.add_argument()
    #Optional arguments
    _opts = _parser.add_argument_group('optional arguments')
    _opts.add_argument()

    #GENENAMECONVERT subparser

    #NORMALIZE subparser

    #SERNASEQ subparser
    se_parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    se_reqs = se_parser.add_argument_group('required arguments')
    se_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory',
        required=True)
    se_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    se_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        required=True)
    #Optional arguments
    se_opts = se_parser.add_argument_group('optional arguments')
    se_opts.add_argument(
        '--output_bed',
        help='Include option to output BED files for each aligned file',
        action='store_true',
        required=False)
    se_opts.add_argument(
        '--output_bigwig',
        help='Include flag to output bigwig files for each aligned file',
        action='store_true',
        required=False)
    se_opts.add_argument(
        '--count_coding',
        help='Include flag to only count reads aligning to protein coding transcripts',
        action='store_true',
        required=False)
    se_opts.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=False)

    #PERNASEQ subparser
    pe_parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    pe_reqs = pe_parser.add_argument_group('required arguments')
    pe_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory -- if paired-end, file names should be exactly the same except for r1/r2.fastq or similar suffix',
        required=True)
    pe_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    pe_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        required=True)
    #Optional arguments
    pe_opts = pe_parser.add_argument_group('optional arguments')
    pe_opts.add_argument(
        '--output_bed',
        help='Include option to output BED files for each aligned file',
        action='store_true',
        required=False)
    pe_opts.add_argument(
        '--output_bigwig',
        help='Include flag to output bigwig files for each aligned file',
        action='store_true',
        required=False)
    pe_opts.add_argument(
        '--count_coding',
        help='Include flag to only count reads aligning to protein coding transcripts',
        action='store_true',
        required=False)
    pe_opts.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=False)

    #RIBOPROF subparser
    rp_parser = subparser.add_parser('', description='', add_help=False)
    #Required arguments
    rp_reqs = rp_parser.add_argument_group('required arguments')
    rp_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory',
        required=True)
    rp_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    rp_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        required=True)
    #Optional arguments
    rp_opts = rp_parser.add_argument_group('optional arguments')
    rp_opts.add_argument(
        '--output_bed',
        help='Include option to output BED files for each aligned file',
        action='store_true',
        required=False)
    rp_opts.add_argument(
        '--output_bigwig',
        help='Include flag to output bigwig files for each aligned file',
        action='store_true',
        required=False)
    rp_opts.add_argument(
        '--count_coding',
        help='Include flag to only count reads aligning to protein coding transcripts',
        action='store_true',
        required=False)
    rp_opts.add_argument(
        '--truncate',
        help='Count reads to truncated reference (truncated reference should be added to the \"reference\" parent directory)',
        action='store_true',
        required=False)
    rp_opts.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=False)

    """
    COLLECT PARSED ARGUMENTS AND PREPARE FOR DOWNSTREAM USE
    """
    #Print help if no arguments/submodules specified
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

    #Parse arguments into NameSpace
    args = parser.parse_args(args)

    #Collect subargs and package, add XPRESSpipe script path to argument dictionary
    args_dict = vars(args)
    args_dict['path'] = __path__

    #Check inputs validity
    check_arguments(args_dict)

    return args, args_dict
