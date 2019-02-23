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

    """
    TRIM SUBPARSER
    """
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
        help='Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if none are provided, software will attempt to auto-detect adaptors -- if \"POLYX\" is provided as a single string in the list, polyX adaptors will be trimmed. If you want to auto-detect adaptors in for paired-end reads, provide a list as ["autoPE"]',
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

    """
    ALIGN SUBPARSER
    """
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
        '-h', '--help',
        action='help',
        help='Show help message and exit')
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
    align_opts.add_argument(
        '-m', '--max_processors',
        help='Number of max processors to use for tasks (default: No limit)',
        metavar='<int>',
        default=DEFAULT_MAX_PROCESSORS,
        required=False)

    """
    COUNT SUBPARSER
    """
    count_parser = subparser.add_parser('count', description='Get counts and a counts table from aligned output', add_help=False)
    #Required arguments
    count_reqs = count_parser.add_argument_group('required arguments')
    count_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory',
        required=True)
    count_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    count_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        required=True)
    #Optional arguments
    count_opts = count_parser.add_argument_group('optional arguments')
    count_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    count_opts.add_argument(
        '--count_coding',
        help='Count reads to coding transcripts only (coding reference should be added to the \"reference\" parent directory)',
        action='store_true',
        required=False)
    count_opts.add_argument(
        '--truncate',
        help='Count reads to truncated reference (truncated reference should be added to the \"reference\" parent directory)',
        action='store_true',
        required=False)
    count_opts.add_argument(
        '--experiment',
        help='Experiment name to name counts table -- if none provided, datetime is used to format file output name',
        type=<str>,
        required=False)

    """
    METAGENE SUBPARSER
    """
    metagene_parser = subparser.add_parser('metagene', description='Compile summarized metagene profiles for each sample in a directory', add_help=False)
    #Required arguments
    metagene_reqs = metagene_parser.add_argument_group('required arguments')
    metagene_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory of SAM alignment files',
        required=True)
    metagene_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    metagene_reqs.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=True)
    metagene_reqs.add_argument(
        '-r', '--reference',
        help='Path to parent organism reference directory',
        type=<str>,
        required=True)
    metagene_reqs.add_argument(
        '-t', '--reference_type',
        help='RefFlat type (i.e. \"DEFAULT\", \"CODING\", \"CODING_TRUNCATED\",)',
        type=<str>,
        required=True)
    #Optional arguments
    metagene_opts = metagene_parser.add_argument_group('optional arguments')
    metagene_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')

    """
    READ DISTRIBUTION SUBPARSER
    """
    distribution_parser = subparser.add_parser('readDistribution', description='Compile summarized distributions for each sample in a directory', add_help=False)
    #Required arguments
    distribution_reqs = distribution_parser.add_argument_group('required arguments')
    distribution_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory of trimmed fastq (or untrimmed fastq) files',
        required=True)
    distribution_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    distribution_reqs.add_argument(
        '--experiment',
        help='Experiment name to name counts table -- if none provided, datetime is used to format file output name',
        type=<str>,
        required=True)
    #Optional arguments
    distribution_opts = distribution_parser.add_argument_group('optional arguments')
    distribution_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')

    """
    PERIODICITY SUBPARSER
    """
    period_parser = subparser.add_parser('periodicity', description='Calculate periodicity of transcripts using the most abundant transcript length for alignments to map per sample', add_help=False)
    #Required arguments
    period_reqs = period_parser.add_argument_group('required arguments')
    period_reqs.add_argument()
    period_reqs.add_argument(
        '-i', '--input',
        help='Path to input directory of SAM alignment files',
        required=True)
    period_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory',
        required=True)
    period_reqs.add_argument(
        '-g', '--gtf',
        help='Path and file name of reference GTF',
        type=<str>,
        required=True)
    period_reqs.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=True)
    #Optional arguments
    period_opts = period_parser.add_argument_group('optional arguments')
    period_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    period_opts.add_argument(
        '--landmark',
        help='Metagene count landmark variable for landmark anchor of periodicity analysis (default: %s)' % 'cds_start',
        default='cds_start',
        metavar='<str>',
        required=False)
    period_opts.add_argument(
        '--downstream',
        help='Number of nucleotides to track after the landmark (default: %s)' % 200,
        default=200,
        metavar='<int>',
        required=False)

    """
    TRUNCATE SUBPARSER
    """
    truncate_parser = subparser.add_parser('truncate', description='Create a coding-only and coding-only truncated reference GTF', add_help=False)
    #Required arguments
    truncate_reqs = truncate_parser.add_argument_group('required arguments')
    truncate_reqs.add_argument(
        '-g', '--gtf',
        help='Path and file name of GTF reference file to process',
        type=<str>,
        required=True)
    #Optional arguments
    truncate_opts = truncate_parser.add_argument_group('optional arguments')
    truncate_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    truncate_opts.add_argument(
        '-t', '--truncate_amount',
        help='Path and file name of GTF reference file to process (default: %s)' % 45,
        default=45,
        type=<int>,
        required=False)
    truncate_opts.add_argument(
        '-c', '--create_refFlats',
        help='Provide flag to output refFlat files for each transcript reference created',
        action='store_true',
        required=False)

    """
    MAKEFLAT SUBPARSER
    """
    makeflat_parser = subparser.add_parser('makeFlat', description='Create flattened reference GTFs', add_help=False)
    #Required arguments
    makeflat_reqs = makeflat_parser.add_argument_group('required arguments')
    makeflat_reqs.add_argument(
        '-i', '--input',
        help='Path where input transcripts*.gtf files are found',
        type=<str>,
        required=True)

    """
    REFERENCE SUBPARSER
    """
    reference_parser = subparser.add_parser('createReference', description='Create a STAR reference directory', add_help=False)
    #Required arguments
    reference_reqs = reference_parser.add_argument_group('required arguments')
    reference_reqs.add_argument(
        '-o', '--output',
        help='Path to directory for output',
        type=<str>,
        required=True)
    reference_reqs.add_argument(
        '-f', '--fasta',
        help='Path to directory containing genomic fasta files',
        type=<str>,
        required=True)
    reference_reqs.add_argument(
        '-g', '--gtf',
        help='Path and file name of reference GTF',
        type=<str>,
        required=True)
    #Optional arguments
    reference_opts = reference_parser.add_argument_group('optional arguments')
    reference_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    reference_opts.add_argument(
        '-t', '--threads',
        help='Specify number of threads to use (default: %s)' % 8,
        default=8,
        type=<int>,
        required=False)
    reference_opts.add_argument(
        '--sjdbOverhang',
        help='Specify number of threads to use (default: %s)' % 100,
        default=100,
        type=<int>,
        required=False)

    """
    RRNAPROBE SUBPARSER
    """
    probe_parser = subparser.add_parser('rrnaProbe', description='Collect overrepresented sequences from multiple FastQC zipfile outputs (IMPORTANT: Run a BLAST search on sequences you choose to use as depletion probes to verify they are rRNAs)', add_help=False)
    #Required arguments
    probe_reqs = probe_parser.add_argument_group('required arguments')
    probe_reqs.add_argument(
        "-i", "--input",
        help="Path to zipped files ",
        metavar='<string>',
        required=True)
    probe_reqs.add_argument(
        "-o", "--output",
        help="Output file name to write output to",
        metavar='<string>',
        required=True)
    #Optional arguments
    probe_opts = probe_parser.add_argument_group('optional arguments')
    probe_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    probe_opts.add_argument(
        '-m', '--min_overlap',
        help='Minimum number of bases that must match on a side to combine sequences (default: %s)' % 5,
        default=5,
        metavar='<int>',
        required=False)
    probe_opts.add_argument(
        '--footprint_only',
        help='Only take zip files that are ribosome profiling footprints (file names must contain \"FP\", \"RPF\", or \"FOOTPRINT\")',
        action='store_true',
        required=False)

    """
    CONVERTNAMES SUBPARSER
    """
    convert_parser = subparser.add_parser('convertNames', description='Convert gene identifiers using a GTF reference', add_help=False)
    #Required arguments
    convert_reqs = convert_parser.add_argument_group('required arguments')
    convert_reqs.add_argument(
        '-d', '--data',
        help='Path and file name to sequence dataframe',
        metavar='<string>',
        required=True)
    convert_reqs.add_argument(
        '-g', '--gtf',
        help='Path and file name to reference GTF',
        metavar='<string>',
        required=True)
    #Optional arguments
    convert_opts = convert_parser.add_argument_group('optional arguments')
    convert_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    convert_opts.add_argument(
        '--orig_name_label',
        help='Label of original name (usually \"gene_id \")',
        default='gene_id \"',
        metavar='<str>',
        required=False)
    convert_opts.add_argument(
        '--orig_name_location',
        help='Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)',
        default=0,
        metavar='<int>',
        required=False)
    convert_opts.add_argument(
        '--new_name_label',
        help='Label of original name (usually \"gene_name \")',
        default='gene_name \"',
        metavar='<str>',
        required=False)
    convert_opts.add_argument(
        '--new_name_location',
        help='Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)',
        default=2,
        metavar='<int>',
        required=False)
    convert_opts.add_argument(
        '--refill',
        help='In some cases, where common gene names are unavailable, the dataframe will fill the gene name with the improper field of the GTF. In this case, specify this improper string and these values will be replaced with the original name',
        default=None,
        metavar='<str>',
        required=False)

    """
    NORMALIZE SUBPARSER
    """
    normalize_parser = subparser.add_parser('normalizeMatrix', description='Perform normalization on sequence matrix', add_help=False)
    #Required arguments
    normalize_reqs = normalize_parser.add_argument_group('required arguments')
    normalize_reqs.add_argument(
        '-d', '--data',
        help='Path and file name to sequence dataframe',
        metavar='<string>',
        required=True)
    #Optional arguments
    normalize_opts = normalize_parser.add_argument_group('optional arguments')
    normalize_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    normalize_opts.add_argument(
        '-m', '--method',
        help='Normalization method to perform (options: \"RPM\", \"RPKM\", \"FPKM\", \"LOG\") -- if using either RPKM or FPKM, a GTF reference file must be included',
        metavar='<string>',
        required=False)
    normalize_opts.add_argument(
        '-g', '--gtf',
        help='Path and file name to reference GTF',
        metavar='<string>',
        required=False)
    normalize_opts.add_argument(
        '--batch',
        help='Include path and filename of dataframe with batch normalization parameters',
        metavar='<string>',
        required=False)

    """
    SERNASEQ SUBPARSER
    """
    se_parser = subparser.add_parser('seRNAseq', description='Trim, align, count, and perform quality control on single-end RNAseq raw data', add_help=False)
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
    se_reqs.add_argument(
        '-e', '--experiment',
        help='Experiment name',
        type=<str>,
        required=True)
    #Optional arguments
    se_opts = se_parser.add_argument_group('optional arguments')
    se_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    se_opts.add_argument(
        '-a', '--adaptors',
        metavar='<list>',
        nargs='+',
        default=None,
        help='Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if none are provided, software will attempt to auto-detect adaptors -- if \"POLYX\" is provided as a single string in the list, polyX adaptors will be trimmed',
        required=False)
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
        '-m', '--method',
        help='Normalization method to perform (options: \"RPM\", \"RPKM\", \"FPKM\", \"LOG\")',
        metavar='<string>',
        required=False)
    se_opts.add_argument(
        '--batch',
        help='Include path and filename of dataframe with batch normalization parameters',
        metavar='<string>',
        required=False)

    """
    PERNASEQ SUBPARSER
    """
    pe_parser = subparser.add_parser('peRNAseq', description='Trim, align, count, and perform quality control on paired-end RNAseq raw data', add_help=False)
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
    pe_reqs.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=True)
    #Optional arguments
    pe_opts = pe_parser.add_argument_group('optional arguments')
    pe_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    pe_opts.add_argument(
        '-a', '--adaptors',
        metavar='<list>',
        nargs='+',
        default=None,
        help='Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if none are provided, software will attempt to auto-detect adaptors -- if \"POLYX\" is provided as a single string in the list, polyX adaptors will be trimmed. If you want to auto-detect adaptors in for paired-end reads, provide a list as ["autoPE"]',
        required=False)
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
        '-m', '--method',
        help='Normalization method to perform (options: \"RPM\", \"RPKM\", \"FPKM\", \"LOG\")',
        metavar='<string>',
        required=False)
    pe_opts.add_argument(
        '--batch',
        help='Include path and filename of dataframe with batch normalization parameters',
        metavar='<string>',
        required=False)

    """
    RIBOPROF SUBPARSER
    """
    rp_parser = subparser.add_parser('riboprof', description='Trim, align, count, and perform quality control on single-end Ribosome Profiling raw data', add_help=False)
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
    rp_reqs.add_argument(
        '--experiment',
        help='Experiment name',
        type=<str>,
        required=True)
    #Optional arguments
    rp_opts = rp_parser.add_argument_group('optional arguments')
    rp_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    rp_opts.add_argument(
        '-a', '--adaptors',
        metavar='<list>',
        nargs='+',
        default=None,
        help='Specify adaptor(s) in list of strings -- if more than one is provided, it will be assumed reads are paired-end -- if none are provided, software will attempt to auto-detect adaptors -- if \"POLYX\" is provided as a single string in the list, polyX adaptors will be trimmed',
        required=False)
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
        '-m', '--method',
        help='Normalization method to perform (options: \"RPM\", \"RPKM\", \"FPKM\", \"LOG\")',
        metavar='<string>',
        required=False)
    rp_opts.add_argument(
        '--batch',
        help='Include path and filename of dataframe with batch normalization parameters',
        metavar='<string>',
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
