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
import time


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

def msg_complete():
    print('Process complete.')

#trim submodule
def msg_trim():
    print("""\n************************\nXPRESSpipe initialized...
        \nAdaptor and quality trimming...\n
  _       ,/'
 (_).  ,/'
  __  ::
 (__)'  `\.
           `\.
        \n************************\n""")

    time.sleep(2)

#align submodule
def msg_align():
    print("""\n************************\nTrimming complete
        \nAligning...\n
                                                    __________
|             _,---.                          _..--'/          \              |
|          ,-'      )                    _,-,'     /            \             |
|         (          )               _,-'  /      /              \            |
|          `-.__, -'             _,-/     /     _/ ,------------. \           |
|     ,--._,---.             _,-'  /     /_,--'' | |     Who    | |           |
|   _(_,-'      )        _,-'/    /_, - '|       | |     is     | |           |
| ,'  (        )      ,'/   /_, -'|      ||_  -  | |    John    | |           |
|(     `--._,-'    _,' /_,-'|     ||_  - |       | |    Galt    | |           |
| `--.__,-'     _,'/,-'|  - ||_ - |      ||_  -  | |      ?     | |           |
|            _,'/,'|  - ___ |     ||_  - |       | `------------' |           |
|        _ ,' |_' -|- _(   ) |_   |      ||_  -  |                |           |
|       ( ): ((`) -| (. `-/ )   - ||_  - |       |                |           |
|      ` | '(-.,') '( _\`/-') _ - |      | _   - |                |           |
|        `   (_.) -(_._ \|,-_)    ||_  - |       |                |           |
| .        `. || ` | (_\||_)  _ - |      ||_  -  |                |           |
|             `.    ` . ||  |   - ||_  - |       |                |           |
|                `.     ||-.|  :  |  _   ||_   - |     __________ |           |
|                   `.  ||    ` -.|   |  |   _   |    | _Doener_ |`'--.._     |
|                      `.           ` - .|  | |  |    || ()._o  ||   __ |     |
|\           `            `.   `           ` -' .|____|`'---.|>_||  |. ||     |
|                            `.                         `' -- .._|__|__||     |
|                      ,----.---------- .                            .-_.     |
|   _                 |.`.,' `.           `.                         || |     |
|    \                \ _|`. / `.            `.                       |_o     |
|     \                | \  `.   \.-----_-----.\     `                        |
|      \               `.|.`| `./ \\  c__)___/ \ \\             .               |
|       \        `         `.   `. \\____\___)__\ \\                            |
|                            `. _ `.\`____________\                           |
|                              | \  |_|   ____   |_|                          |
|          _                   `.|\ |#|__[jrei]__|_|`.                        |
|           \                      `================'  `.                     |


        \n************************\n""")

    time.sleep(2)

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
    time.sleep(2)

#Normalize
def msg_normalize():
    print("""\n************************\nCount tables generated
        \nNormalizing counts based on parameters provided...\n
                           ___,
                    o___.-' /
                    |      _\_
                    |___.-'   `
                    |
                    |
            _   _   j   _   _
           [_]_[_]_[_]_[_]_[_]
           [__j__j__j__j__j__]
             [_j__j__j__j__]
             [__j__j__j__j_]
             [_j__j/V\_j__j]
             [__j_// \\__j_]
             [_j__|   |_j__]
             [__j_|___|__j_]
             [_j__j__j__j__]
             [__j__j__j__j_]
  _   _   _  [_j__j__j__j__]  _   _   _   _
_[_]_[_]_[_]_[__j__j__j__j_]_[_]_[_]_[_]_[_]_
  _j__j__j__j[_j__j__j__j__]j__j__j__j__j_
     j  j  j [  j  j  j  j ] j  j  j  j
        \n************************\n""")

    time.sleep(2)

#Quality control
def msg_quality():
    print("""
    \n****************\nNormalization complete\n****************\n
    \n*********************\nPerforming quality control on the processed sequencing data...\n*********************\n
      _               _               _
     | |             | |             | |
  ___| |__   ___  ___| | _____  _   _| |_
 / __| '_ \ / _ \/ __| |/ / _ \| | | | __|
| (__| | | |  __/ (__|   < (_) | |_| | |_
 \___|_| |_|\___|\___|_|\_\___/ \__,_|\__|
 \n***********************************************************************\n""")

    time.sleep(2)

#final message
def msg_finish():
    print("""
    \n************************\nQuality control complete
    \n************************\nXPRESSpipe processing complete
    \n*******************\n
        """)
