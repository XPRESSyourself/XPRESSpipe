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

__version__ = '0.1.4-beta'

from .__main__ import *
from .align import *
from .arguments import *
from .buildCommand import *
from .compile import *
from .complexity import *
from .convert import *
from .count import *
from .geneCoverage import *
from .gtfModify import *
from .gtfTruncate import *
from .messages import *
from .metagene import *
from .normalizeMatrix import *
from .parallel import *
from .periodicity import *
from .processBAM import *
from .quality import *
from .readDistribution import *
from .rrnaProbe import *
from .test import *
from .trim import *
from .utils import *
