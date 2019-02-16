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

#All ASCII art from http://ascii.co.uk/art

"""
FUNCTIONS
"""
#show license information
def msg_license():
    print("""
RiboPipe
An assembly and analysis pipeline for sequencing data
alias: ripopipe

Copyright (C) 2018  Jordan A. Berg
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
    """)

#trim submodule
def msg_trim():
    print("""\n************************\nRiboPipe initialized...
        \nAdaptor and quality trimming...\n
  _       ,/'
 (_).  ,/'
  __  ::
 (__)'  `\.
           `\.
        \n************************\n""")

#align submodule
def msg_align():
    print("""\n************************\nTrimming complete
        \nAligning...\n
        \n************************\n""")

#count submodule
def msg_count():
    print("""\n************************\nAlignment complete
 ______________________________________________________________________________
|==============================================================================|
|NOVA TOTIUS TERRARUM ORBIS GEOGRAPHICA AC HYDROGRAPHICA TABULA auct:Iud:Hondio|
|===,%% =========|====_== ,%% ================================================_|
|   %%'}  VER  -/+\- (-(  ''%% TERRA _,---._ AQUA~ O ~~~~~~~~   AESTAS    %% /_|
|  ,-' '-.  .oO |_|   `-)_/\..`-_  ,''._,-.)`.~~~_/ /> ~~~~~    _,--.(>  %''%|_|
| /( \(.) \_`|'_,---'''''---`-.-.\/-.   ,  o,-\,',-',---'''''---.`,,' o ,') (> |
| \\) `( `--,-'              (`-. |-=o (( (|)-| ,-;, ,_`'-._/|_   `-. |) /   \  |
|  /    \ ,'              ,_ _)O `\ ,-".));`-./' /_;_|\`,-,    \  _\_`ll |___| |
|  |____,' \             (| `._    `._. \ ._,' _,--'  `'''      \|  \`-`. ___  |
| _ |/|/    `-.____            )     \`---'/  / ,--..`.,-,  _       _\  _\\_/ \ |
|(+)| /            ( \   .---,;. MAR  \   / _,`-.._'`' ,\  `-'      |\|;_|\.__,|
| _,-/              \|\  |  ..--. DEL :\ / /       '--../ `. _       \ ;   \`--|
|/ ,:                `\'-.`_-_ -.. NORT : |             \\,'  `. ,'. | :   .:    |
|.| |        MAR DEL       ,' '-._      |  \____      O ,'     |/  |-' ;   :|  |
| ,-|---------------------|-------',__--|-------\------/-------'--\'-()--.--|`. |
|' ,|            .. ZUR   \           \ |OCEANUS|  , o |    MAR DI `_ /'-   |\ |
| / :    ` ..              |          / : AETHIO \     | |   INDIA     `'   :` |
|:|  \                     |       _,' / \ PICUS  |   /  '        _        /.\:|
|AUTUM\                    |   _;-'   /   \       \_,'           / \_  /  /-.| |
|o|NUS \   TERRA         ._|,-'      /,---.\                  _.'    `'  / `.`o|
| | \_O_`.   AUSTRALIS     `._     ,' ,-. ) `.   __......____/TERRA    ,'HYEMS |
|.\ _|_/> `.   INCOGNITA      `._,/ -'   /_\  \./        AUSTRALIS   ,'    0 |:|
|   \ /     _-._             _,-' |,-.-.    o | `-._   INCOGNITA _,-'  )  /_\|.|
|--. / . : /____\---.....---',\_O_\   ,-. -')>/\_O_ `---.....---'('')  (  /_/  |
|--'  `.\;  |--|-. ____IGNIS  _|_/|`._>' `-_,'   \)`-   AER `-' ,  (""\  \     |
|______________________________________________________________________________|
        \nGenerating count tables...\n
        1...2...3...4...5...6...7...8...9...
        \n************************\n""")

#check output submodule
def msg_checking():
    print("""\n************************\nCount tables generated
        \nGenerating quality control summaries...\n
      _               _               _
     | |             | |             | |
  ___| |__   ___  ___| | _____  _   _| |_
 / __| '_ \ / _ \/ __| |/ / _ \| | | | __|
| (__| | | |  __/ (__|   < (_) | |_| | |_
 \___|_| |_|\___|\___|_|\_\___/ \__,_|\__|
        \n************************\n""")

#clean submodule
def msg_cleaning():
    print("""
    \n****************\nQuality control summaries complete\n****************\n
    \n*********************\nCleaning up the output...\n*********************\n
    .-.
    | |
    |=|
    |=|
    | |
    | |
    | |
    | |
    | |
    | |
    | |
    | |
    | |
    | |
    | |
    |=|
    |=|
    |_|
  .=/I\=.
 ////V\\\\
 |#######|
 |||||||||
 |||||||||
 |||||||||
 \n***********************************************************************\n""")

#final message
def msg_finish():
    print("""
    \n************************\nRiboPipe processing complete
    \nIt is recommended that you check the number of uniquely mapped reads in
    each of the samples. Samples with less than 2 million mapped reads are often
    considered to not allow for accuracy in downstream analyses. This number can
    vary from organism to organism.
    This information can be found in the
    .../outputDir/assembly/counts/XXX_raw_counts_compiled.csv output file.
    \nPlease refer to the diffex sub-module for differential expression analysis: 'ribopipe diffex --help'
    \nOther programs for differential expression analysis in ribosome profiling
    data include Xtail and RiboDiff
    https://www.ncbi.nlm.nih.gov/pubmed/27041671
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5198522/\n*******************\n
        """)
