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
from .__init__ import __version__
from .messages import *
from .arguments import get_arguments
from .trim import run_trim
from .align import run_align, run_seRNAseq, run_peRNAseq, run_riboprof
from .count import run_count
from .quality import run_quality, read_distribution, metagene, periodicity
from .reference import run_truncate, run_rrnaprobe

"""
DESCRIPTION: Main function to call necessary functions for sub-modules

ASSUMPTIONS:
Proper arguments are provided where some user renaming of files may be required
"""
def main(args=None):

    #Print license information
    msg_license()

    #Collate CLI arguments provided by user
    try:
        args, args_dict = get_arguments(args, __version__)
    except:
        raise Exception("There was an issue in processing the arguments for the provided function.")

    #Execute corresponding functions determined by arguments provided by user
    if args.cmd == 'trim':
        run_riboprof(args_dict)

    elif args.cmd == 'align':
        run_riboprof(args_dict)

    elif args.cmd == 'count':
        run_riboprof(args_dict)

    elif args.cmd == 'quality':
        run_riboprof(args_dict)

    elif args.cmd == 'truncate':
        run_riboprof(args_dict)

    elif args.cmd == 'rrnaprobe':
        run_riboprof(args_dict)

    elif args.cmd == 'seRNAseq':
        run_riboprof(args_dict)

    elif args.cmd == 'peRNAseq':
        run_riboprof(args_dict)

    elif args.cmd == 'riboprof':
        run_riboprof(args_dict)

    else:
        raise Exception("Invalid function processing function provided.")

"""
DESCRIPTION: Run main
"""
if __name__ == "__main__":

    sys.exit(main() or 0)
