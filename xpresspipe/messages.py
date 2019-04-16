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
"""IMPORT DEPENDENCIES"""
import os
import sys
import time

"""Check progress report for errors and exceptions"""
def check_process(
        log_file,
        message_func,
        step):

    try:
        os.system(
            '[[ $(cat ' + log_file + ' | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]]'
            + ' || { echo "Errors or exceptions were present in ' + step + ', please refer to the '
            + str(log_file[log_file.rfind('/') + 1:])
            + ' file for information concerning errors"; exit 1; }')
    except:
        message_func

"""Show license information"""
def msg_license():
    print("""
XPRESSpipe
An assembly and analysis pipeline for sequencing data
alias: xpresspipe

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

"""Print completed process message"""
def msg_complete():
    print("""Process complete.
************************
    """)

"""Trim submodule message"""
def msg_trim():
    print("""\nXPRESSpipe initialized...
    \nAdaptor and quality trimming...\n
  _       ,/'
 (_).  ,/'
  __  ::
 (__)'  `\.
           `\.
       """)

    time.sleep(1)

"""Align submodule message"""
def msg_align():
    print("""Aligning...\n
                                                    __________
|             _,---.                          _..--'/                       |
|          ,-'      )                    _,-,'     /                        |
|         (          )               _,-'  /      /                         |
|          `-.__, -'             _,-/     /     _/ ,------------.           |
|     ,--._,---.             _,-'  /     /_,--'' | |            | |           |
|   _(_,-'      )        _,-'/    /_, - '|       | |   XPRESS   | |           |
| ,'  (        )      ,'/   /_, -'|      ||_  -  | |            | |           |
|(     `--._,-'    _,' /_,-'|     ||_  - |       | |  yourself  | |           |
| `--.__,-'     _,'/,-'|  - ||_ - |      ||_  -  | |            | |           |
|            _,'/,'|  - ___ |     ||_  - |       | `------------' |           |
|        _ ,' |_' -|- _(   ) |_   |      ||_  -  |                |           |
|       ( ): ((`) -| (. `-/ )   - ||_  - |       |                |           |
|      ` | '(-.,') '( _\`/-') _ - |      | _   - |                |           |
|        `   (_.) -(_._|,-_)    ||_  - |       |                |           |
| .        `. || ` | (_\||_)  _ - |      ||_  -  |                |           |
|             `.    ` . ||  |   - ||_  - |       |                |           |
|                `.     ||-.|  :  |  _   ||_   - |     __________ |           |
|                   `.  ||    ` -.|   |  |   _   |    | ________ |`'--.._     |
|                      `.           ` - .|  | |  |    || ()._o  ||   __ |     |
|\           `            `.   `           ` -' .|____|`'---.|>_||  |. ||     |
|                            `.                         `' -- .._|__|__||     |
|                      ,----.---------- .                            .-_.     |
|   _                 |.`.,' `.           `.                         || |     |
|                   _|`. / `.            `.                       |_o     |
|                    |  `.  .-----_-----.\     `                        |
|                    `.|.`| `./\  c__)___/\             .               |
|              `         `.   `.\____\___)__\\                            |
|                            `. _ `.\`____________\                           |
|                              |  |_|   ____   |_|                          |
|          _                   `.|\ |#|__[....]__|#|`.                        |
|                                `================'  `.                     |


       """)

    time.sleep(1)

"""Count submodule message"""
def msg_count():
    print("""Generating count tables...\n
    1...2...3...4...5...6...7...8...9...
 ______________________________________________________________________________
|==============================================================================|
|NOVA TOTIUS TERRARUM ORBIS GEOGRAPHICA AC HYDROGRAPHICA TABULA auct:Iud:Hondio|
|===,%% =========|====_== ,%% ================================================_|
|   %%'}  VER  -/+\- (-(  ''%% TERRA _,---._ AQUA~ O ~~~~~~~~   AESTAS    %% /_|
|  ,-' '-.  .oO |_|   `-)_/\..`-_  ,''._,-.)`.~~~_/ /> ~~~~~    _,--.(>  %''%|_|
| /((.)_`|'_,---'''''---`-.-.\/-.   ,  o,-\,',-',---'''''---.`,,' o ,') (> |
|\) `( `--,-'              (`-. |-=o (( (|)-| ,-;, ,_`'-._/|_   `-. |) /    |
|  /    ,'              ,_ _)O `\ ,-".));`-./' /_;_|\`,-,     _\_`ll |___| |
|  |____,'             (| `._    `._. ._,' _,--'  `'''     | `-`. ___  |
| _ |/|/    `-.____            )    `---'/  / ,--..`.,-,  _       _\  _\\_/ |
|(+)| /            (   .---,;. MAR    / _,`-.._'`' ,\  `-'      |\|;_|\.__,|
| _,-/             |\  |  ..--. DEL :\ / /       '--../ `. _       ;  `--|
|/ ,:                `\'-.`_-_ -.. NORT : |            \,'  `. ,'. | :   .:    |
|.| |        MAR DEL       ,' '-._      | ____      O ,'     |/  |-' ;   :|  |
| ,-|---------------------|-------',__--|-------\------/-------'--\'-()--.--|`. |
|' ,|            .. ZUR             |OCEANUS|  , o |    MAR DI `_ /'-   |\ |
| / :    ` ..              |          / : AETHIO     | |   INDIA     `'   :` |
|:|                      |       _,' / PICUS  |   /  '        _        /.\:|
|AUTUM\                    |   _;-'   /        _,'           /_  /  /-.| |
|o|NUS   TERRA         ._|,-'      /,---.\                  _.'    `'  / `.`o|
| |_O_`.   AUSTRALIS     `._     ,' ,-. ) `.   __......____/TERRA    ,'HYEMS |
|.\ _|_/> `.   INCOGNITA      `._,/ -'   /_\ ./        AUSTRALIS   ,'    0 |:|
|   /     _-._             _,-' |,-.-.    o | `-._   INCOGNITA _,-'  )  /_\|.|
|--. / . : /____\---.....---',\_O_\   ,-. -')>/\_O_ `---.....---'('')  (  /_/  |
|--'  `.\;  |--|-. ____IGNIS  _|_/|`._>' `-_,'  )`-   AER `-' ,  (""\      |
|______________________________________________________________________________|

       """)
    time.sleep(1)

"""Normalize submodule message"""
def msg_normalize():
    print("""Normalizing counts based on parameters provided...\n
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
             [__j_//\__j_]
             [_j__|   |_j__]
             [__j_|___|__j_]
             [_j__j__j__j__]
             [__j__j__j__j_]
  _   _   _  [_j__j__j__j__]  _   _   _   _
_[_]_[_]_[_]_[__j__j__j__j_]_[_]_[_]_[_]_[_]_
  _j__j__j__j[_j__j__j__j__]j__j__j__j__j_
     j  j  j [  j  j  j  j ] j  j  j  j
       """)

    time.sleep(1)

"""Quality control submodule message"""
def msg_quality():
    print("""Performing quality control on the processed sequencing data...\n
      _               _               _
     | |             | |             | |
  ___| |__   ___  ___| | _____  _   _| |_
 / __| '_ / _/ __| |/ / _| | | | __|
| (__| | | |  __/ (__|   < (_) | |_| | |_
___|_| |_|\___|\___|_|\_\___/__,_|\__|
\n""")

    time.sleep(1)

"""Final message"""
def msg_finish():
    print("""Quality control complete\n
             __
            / /
           / /
          / /
 __      / /
 \ \    / /
  \ \  / /
   \ \/ /
    \__/

   XPRESSpipe processing complete

                                    ..
                                     .(  )`-._
                                   .'  ||     `._
                                 .'    ||        `.
                              .'       ||          `._
                            .'        _||_            `-.
                         .'          |====|              `..
                       .'            __/               (  )
                     ( )               ||          _      ||
                     /|\               ||       .-`     ||
                   .' | '              ||   _.-'    |     ||
                  /   |\             || .'   `.__.'     ||   _.-..
                .'   /| `.            _.-'   _.-'       _.-.`-'`._`.`
                 .' |  |        .-.`    `./      _.-`.    `._.-'
                 |.   |  `.   _.-'   `.   .'     .'  `._.`---`
                .'    |   |  :   `._..-'.'        `._..'  ||
               /      |    `-._.'    ||                 ||
              |     .'|`.  |           ||_.--.-._         ||
              '    /  |        __.--'\    `. :        ||
                .'  |  |   ..-'     `._-._.'        ||
`.._            |/    |    `.      `._.-              ||
    `-.._       /     |       `-.'_.--'                 ||
         `-.._.'      |      |        | |         _ _ _  _'_ _ _ _ _
              `-.._   |             | |        |_|_|_'|_|_|_|_|_|_|
                  [`--^-..._.'        | |       /....../|  __   __  |
                  `---.._|`--.._    | |      /....../ | |__| |__| |
                   __  _ `-.._| `-._|_|_ _ _/_ _ _ /  | |__| |__| |
                       _o_   _`-._|_|_|_|_|_|_|_|_/   '-----------/
                     _`.|.'  _  - .--.--.--.--.--.`--------------'
      .```-._ ``-.._  __   _    _ '--'--'--'--'--'  - _ - _  __/
 .`-.```-._ ``-..__``.- `.      _     -  _  _  _ -    _-   _  __/(.``-._
 _.-` ``--..  ..    _.-` ``--..  .. .._ _. __ __ _ __ ..--.._ / .( _..``
`.-._  `._  `- `-._  .`-.```-._ ``-..__``.-  -._--.__---._--..-._`...```
   _.-` ``--..  ..  `.-._  `._  `- `-._ .-_. ._.- -._ --.._`` _.-`---`-.



        """)
