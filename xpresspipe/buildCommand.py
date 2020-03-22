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
from __future__ import print_function

import os
import sys
import numpy as np
from math import ceil, floor

# ASCII art from David Palmer https://www.asciiart.eu/people/faces
def welcome():
    print("""
                _____    ____
             .#########.#######..
          .#######################.
        .############################.
       .################################.
      .#########,##,#####################.
     .#########-#,'########## ############.
    .######'#-##' # ##'### #. `####:#######.
    #####:'# #'###,##' # # +#. `###:':######
   .####,'###,##  ###  # # #`#. -####`######.
  .####+.' ,'#  ##' #   # # #`#`.`#####::####
  ####'    #  '##'  #   #_'# #.## `#######;##
  #:##'      '       #   # ``..__# `#######:#
  #:##'  .#######s.   #.  .s######.\`#####:##
  #:##   ."______""'    '""_____"". `.#####:#
 .#:##   ><'(##)'> )    ( <'(##)`><   `####;#
 ##:##  , , -==-' '.    .` `-==- . \  ######'
 #|-'| / /      ,  :    : ,       \ ` :####:'
 :#  |: '  '   /  .     .  .  `    `  |`####
 #|  :|   /   '  '       `  \   . ,   :  #:#        WELCOME TO THE XPRESSPIPE
 #L. | | ,  /   `.-._ _.-.'   .  \ \  : ) J##     / COMMAND BUILDER! I will
###\ `  /  '                   \  : : |  /##     /  be your guide to help you
 ## #|.:: '       _    _        ` | | |'####    /   through the necessary steps
 #####|\  |  (.-'.__`-'__.`-.)   :| ' ######        to build your reference or
 ######\:      `-...___..-' '     :: /######        to process your sequence data.
 #######\``.                   ,'|  /#######
.# ######\  \       ___       / /' /#######
# #'#####|\  \    -     -    /  ,'|### #. #.
`#  #####| '-.`             ' ,-'  |### #  #
    #' `#|    '._         ,.-'     #`#`#.
         |       .'------' :       |#   #
         |       :         :       |
         |       :         :       |
                 :         :          dp


    """)

# ASCII art from David Palmer https://www.asciiart.eu/people/faces
def ref_guide():
    print("""
                _____    ____
             .#########.#######..
          .#######################.
        .############################.
       .################################.
      .#########,##,#####################.
     .#########-#,'########## ############.
    .######'#-##' # ##'### #. `####:#######.
    #####:'# #'###,##' # # +#. `###:':######
   .####,'###,##  ###  # # #`#. -####`######.
  .####+.' ,'#  ##' #   # # #`#`.`#####::####
  ####'    #  '##'  #   #_'# #.## `#######;##
  #:##'      '       #   # ``..__# `#######:#
  #:##'  .#######s.   #.  .s######.\`#####:##
  #:##   ."______""'    '""_____"". `.#####:#
 .#:##   ><'(##)'> )    ( <'(##)`><   `####;#
 ##:##  , , -==-' '.    .` `-==- . \  ######'
 #|-'| / /      ,  :    : ,       \ ` :####:'
 :#  |: '  '   /  .     .  .  `    `  |`####
 #|  :|   /   '  '       `  \   . ,   :  #:#        Sweet! Let's build a reference!
 #L. | | ,  /   `.-._ _.-.'   .  \ \  : ) J##     / Please be so kind and answer
###\ `  /  '                   \  : : |  /##     /  the following questions for
 ## #|.:: '       _    _        ` | | |'####    /   me to know how to best help you!
 #####|\  |  (.-'.__`-'__.`-.)   :| ' ######
 ######\:      `-...___..-' '     :: /######
 #######\``.                   ,'|  /#######
.# ######\  \       ___       / /' /#######
# #'#####|\  \    -     -    /  ,'|### #. #.
`#  #####| '-.`             ' ,-'  |### #  #
    #' `#|    '._         ,.-'     #`#`#.
         |       .'------' :       |#   #
         |       :         :       |
         |       :         :       |
                 :         :          dp


    """)

# ASCII art from David Palmer https://www.asciiart.eu/people/faces
def run_guide():
    print("""
                _____    ____
             .#########.#######..
          .#######################.
        .############################.
       .################################.
      .#########,##,#####################.
     .#########-#,'########## ############.
    .######'#-##' # ##'### #. `####:#######.
    #####:'# #'###,##' # # +#. `###:':######
   .####,'###,##  ###  # # #`#. -####`######.
  .####+.' ,'#  ##' #   # # #`#`.`#####::####
  ####'    #  '##'  #   #_'# #.## `#######;##
  #:##'      '       #   # ``..__# `#######:#
  #:##'  .#######s.   #.  .s######.\`#####:##
  #:##   ."______""'    '""_____"". `.#####:#
 .#:##   ><'(##)'> )    ( <'(##)`><   `####;#
 ##:##  , , -==-' '.    .` `-==- . \  ######'
 #|-'| / /      ,  :    : ,       \ ` :####:'
 :#  |: '  '   /  .     .  .  `    `  |`####
 #|  :|   /   '  '       `  \   . ,   :  #:#        Awesome! You've chose to
 #L. | | ,  /   `.-._ _.-.'   .  \ \  : ) J##     / build a command to run
###\ `  /  '                   \  : : |  /##     /  your data. Please be so
 ## #|.:: '       _    _        ` | | |'####    /   kind and answer the following
 #####|\  |  (.-'.__`-'__.`-.)   :| ' ######        questions for me to know how
 ######\:      `-...___..-' '     :: /######        to best help you!
 #######\``.                   ,'|  /#######
.# ######\  \       ___       / /' /#######
# #'#####|\  \    -     -    /  ,'|### #. #.
`#  #####| '-.`             ' ,-'  |### #  #
    #' `#|    '._         ,.-'     #`#`#.
         |       .'------' :       |#   #
         |       :         :       |
         |       :         :       |
                 :         :          dp


    """)





def build_curation():

    output = input('Where would you like to output your parent reference directory (provide full path, no spaces)?: ')
    fasta = input('Please provide the full path (no spaces) to the genomic fasta files (ideally, these will be in their own folder): ')
    gtf = input('Please provide the path and filename of the GTF will use: ')

    mod = input('Would you like to modify your GTF for use in quantification during the pipeline? (yes/no): ')
    if mod.lower().replace(' ','') == 'yes':
        mod_value = ''
        add1 = input('Would you like to only quantify to Ensembl canonical transcripts (longest transcript only, not generally necessary)? (yes/no): ')
        if add1.lower().replace(' ','') == 'yes':
            mod_value = mod_value + ' -l'
        add2 = input('Would you like to only quantify to protein coding transcripts? (yes/no): ')
        if add2.lower().replace(' ','') == 'yes':
            mod_value = mod_value + ' -p'
        add3 = input('Are you performing ribosome profiling or would other like to truncate each transcripts CDS records for quantification? (yes/no): ')
        if add3.lower().replace(' ','') == 'yes':
            mod_value = mod_value + ' -t'

        if add3.lower().replace(' ','') == 'yes':
            five_prime = input('Would you like to use a value for 5\' truncation other than the default of 45 nucleotides? (yes/no): ')
            if five_prime.lower().replace(' ','') == 'yes':
                five_prime = input('Please enter the value you would like to use (must be an int): ')
                five_prime = ' --truncate_5prime ' + str(five_prime)
            else:
                five_prime = ''

        three_prime = input('Would you like to use a value for 3\' truncation other than the default of 15 nucleotides? (yes/no): ')
        if three_prime.lower().replace(' ','') == 'yes':
            three_prime = input('Please enter the value you would like to use (must be an int): ')
            three_prime = ' --truncate_3prime ' + str(three_prime)
        else:
            three_prime = ''

    else:
        mod_value = ''
        five_prime = ''
        three_prime = ''

    sjdb = input('What is your read size (is using paired-end sequencing, what is the length of one of the mates)? Input an integer value: ')
    sjdb_value = ' --sjdbOverhang ' + str(int(sjdb) - 1)

    genome = input('Does your organism have a relatively small genome (as compared to Homo sapiens)? (yes/no): ')
    if genome.lower().replace(' ','') == 'yes':
        genome_size = input('Approximately how many bases long is the genome? Input an integer value: ')

        genome_value = ' --genome_size ' + str(genome_size)
    else:
        genome_value = ''

    processors = input('Would you like to throttle number of processors (default will use all available)? (yes/no): ')
    if processors.lower() == 'no':
        proc_value = ''
    else:
        proc_value = input('What is the maximum number of processors you would like to use? (input must be integer): ')
        proc_value = ' --max_processors ' + str(proc_value)

    # Summarize
    print('\nHere is the command you should use to process your data: ')
    print(
        'xpresspipe curateReference'
        + ' -o ' + str(output)
        + ' -f ' + str(fasta)
        + ' -g ' + str(gtf)
        + str(mod_value)
        + str(five_prime)
        + str(three_prime)
        + str(sjdb_value)
        + str(genome_value)
        + str(proc_value)
        )
    print('\nPlease double-check the output for accuracy and that all desired arguments were caught by the builder and incorporated into the command.')

    # Run
    run = input('Do you want to run this command now? (yes/no): ')
    if run.lower().replace(' ','') == 'yes':
        os.system(
            'xpresspipe curateReference'
            + ' -o ' + str(output)
            + ' -f ' + str(fasta)
            + ' -g ' + str(gtf)
            + str(mod_value)
            + str(five_prime)
            + str(three_prime)
            + str(sjdb_value)
            + str(genome_value)
            + str(proc_value)
            )
    else:
        sys.exit(1)


def build_pipeline():

    pipeline = input('What kind of processing would you like to perform (select the number)?\n1: Ribosome Profiling\n2: Single-end RNA-seq\n3: Paired-end RNA-seq\nInput choice: ')

    if pipeline == '1':
        pipeline = 'riboseq'
    elif pipeline == '2':
        pipeline = 'seRNAseq'
    elif pipeline == '3':
        pipeline = 'peRNAseq'
    else:
        print('Invalid input, try again...\n')
        sys.exit(1)

    # Get input directory
    print('Note: There are some arguments that this builder will not include.\nIf you want to use these options, please add them to the command output at the end and run yourself.\n')
    indir = input('Where is your raw sequence data found? (provide full path to directory with .fastq files, no spaces): ')

    # Get output directory
    outdir = input('Where would you like your processed data output to? (provide full path, no spaces): ')

    # References
    reference_dir = input('Where are your genome index and other reference files stored? (provide full path to parent directory): ')
    gtf = input('Which GTF file would you like to use to quantify your alignments with (provide full path and file name, no spaces)?\nExample: /path/to/file/transcripts_LCT.gtf for a modified GTF using all options: ')

    experiment = input('What would you like to name your experiment (no spaces)?: ')

    # Get adaptors
    adaptors_bool = input('Were adaptors used in generating your sequence libraries? (yes/no): ')
    if adaptors_bool.lower().replace(' ','') == 'yes':
        polyx_bool = input('Were polyX adaptors used? (yes/no): ')
        if polyx_bool.lower().replace(' ','') == 'yes':
            if pipeline == 'seRNAseq' or pipeline == 'riboseq':
                adaptors = ' -a polyX'
            else:
                raise Warning('polyX adaptors not compatible with paired-end sequencing module currently. Please try again')
                sys.exit(1)
        else:
            if pipeline == 'seRNAseq' or pipeline == 'riboseq':
                adaptors = input('Please provide the adaptor used: ')
                adaptors = ' -a ' + str(adaptors)
            else:
                adaptors = input('Please provide the adaptors used separated by a space: ')
                if ' ' not in adaptors:
                    raise Warning('Input did not include a space. Please try again')
                    sys.exit(1)
                else:
                    adaptors = ' -a ' + str(adaptors)
    else:
        if pipeline == 'peRNAseq':
            adaptors = ' -a None None'
        else:
            adaptors = ' -a None'

    twopass_bool = input('Would you like to perform two-pass alignment to search for novel splice junctions?\nThis will significantly increase the time it takes to align and is not usually necessary for well-studied organisms (yes/no): ')
    if twopass_bool.lower().replace(' ','') == 'yes':
        twopass = ' --two-pass'
    else:
        twopass = ''

    quality_input = input('Would you like to choose a minimum quality score other than 28? (yes/no) ')
    if str(quality_input).lower().replace(' ','') == 'yes':
        quality_value = input('Please enter the value you would like to use (must be an integer equal to or above 0): ')
        if isinstance(int(quality_value), int) and int(quality_value) >= 0:
            quality = ' -q ' + str(quality_value)
        else:
            quality = ''
            print('Invalid input provided. If you want to include a minimum cut-off, add \"-q your_integer\" to the end of the command the builder gives you.')
    else:
        quality = ''

    length_input = input('Would you like to choose a minimum read length other than 18? (yes/no) ')
    if str(length_input).lower().replace(' ','') == 'yes':
        length_value = input('Please enter the value you would like to use (must be an integer equal to or above 0): ')
        if isinstance(int(length_value), int) and int(length_value) >= 0:
            length = ' -q ' + str(length_value)
        else:
            length = ''
            print('Invalid input provided. If you want to include a minimum cut-off, add \"-l your_integer\" to the end of the command the builder gives you.')
    else:
        length = ''

    front_trim = input('Do you want to perform 5\' trimming of a certain number of bases? (yes/no): ')
    if front_trim.lower().replace(' ','') == 'yes':
        front_value = input('How many bases would you like trimmed? ')
        if isinstance(int(front_value), int) and int(front_value) >= 0:
            front = ' --front_trim ' + str(front_value)
        else:
            front = ''
            print('Invalid input provided. If you want to include 5\' base trimming, add \"--front_trim your_integer\" to the end of the command the builder gives you.')
    else:
        front = ''

    umi_input = input('Did your library preparation method use UMIs? (yes/no): ')
    if umi_input.lower().replace(' ','') == 'yes':
        region = input('Was the UMI in the 5\' or 3\' adapter? (5/3): ')
        if str(region) == '5':
            umi_five_input = input('Please refer to the fastp documentation for UMIs and provide the appropriate option for your usage of UMIs: ')
            fastp_list = [
                'index1',
                'index2',
                'read1',
                'read2',
                'per_index',
                'per_read']

            if umi_five_input in fastp_list:
                umi_loc = ' --umi_location ' + str(umi_five_input)
            else:
                print('Could not find the provided input from fastp\'s list of acceptable inputs.')
                umi_loc = ''
        elif str(region) == '3':
            umi_loc = ' --umi_location 3prime'
        else:
            umi_loc = ''
            print('Invalid input provided. If you want to include a UMI specification, please see the documentation and add the appropriate argument flags and values.')

        if umi_loc != '':
            umi_len_value = input('Please provide the length of your UMI in basepairs: ')
            if isinstance(int(umi_len_value), int) and int(umi_len_value) >= 0:
                umi_len = ' --umi_length ' + str(umi_len_value)
            else:
                umi_len = ''
                print('Invalid input provided. If you want to include a UMI specification, please see the documentation and add the appropriate argument flags and values.')

            spacer_value = input('Please provide the length of your UMI spacer in basepairs. If none was used, provide 0 as your input: ')
            if isinstance(int(spacer_value), int) and int(spacer_value) >= 0:
                umi_spacer = ' --spacer_length ' + str(spacer_value)
            else:
                umi_spacer = ''
                print('Invalid input provided. If you want to include a UMI specification, please see the documentation and add the appropriate argument flags and values.')

        else:
            umi_len = ''
            umi_spacer = ''
    else:
        umi_loc = ''
        umi_len = ''
        umi_spacer = ''

    multimappers_value = input('Would you like to remove multimapped reads before quantification (if a UMI sequence was provided above, this step is unnecessary and you can provide \"no\" as your input)? (yes/no): ')
    if multimappers_value.lower().replace(' ','') == 'yes':
        multimappers = ' --no_multimappers'
    else:
        multimappers = ''

    duplicates_value = input('Would you like to remove duplicated reads before quantification?\nWarning: this may over-compensate and remove biologically relevant reads, especially for shorter reads (if a UMI sequence was provided above, this step is unnecessary and you can provide \"no\" as your input): (yes/no) ')
    if duplicates_value.lower().replace(' ','') == 'yes':
        duplicates = ' --deduplicate'
    else:
        duplicates = ''

    bed_value = input('Would you like to output a BED-formatted file? (yes/no): ')
    if bed_value.lower().replace(' ','') == 'yes':
        bed = ' --output_bed'
    else:
        bed = ''

    remove_rrna_alignments = input('Would you like to remove rRNA alignments from the alignment files? (yes/no): ')
    if remove_rrna_alignments.lower().replace(' ','') == 'yes':
        rrna_alignments = ' --remove_rrna'
    else:
        rrna_alignments = ''

    quant = input('By default, XPRESSpipe uses HTSeq for read quantification.\nHowever, if you would like to estimate isoform abundances, we suggest you use Cufflinks.\nWould you like to use Cufflinks? (yes/no): ')
    if quant.lower().replace(' ','') == 'yes':
        quant_value = ' --quantification_method cufflinks'
    else:
        quant_value = ''

    feature = input('By default, if using HTSeq, quantifications will be made to the exons records, unless working with ribosome profiling data, in which case the CDS will be used.\nWould you like to change this default? (yes/no): ')
    if feature.lower().replace(' ','') == 'yes':
        feature_value = input('Please provide the feature you would like to quantify to, and make sure your entry is case-sensitive: ')
        feature_value = feature_value.replace(' ','')
        feature_value = ' --feature_type ' + str(feature_value)
    else:
        feature_value = ''

    stranded = input('Did you use a strand-aware RNA-seq kit? (yes/no): ')
    if stranded.lower().replace(' ','') == 'yes':
        if quant_value == 'htseq':
            stranded_value = ' --stranded yes'
        else:
            stranded_value = input('Please input the appropriate Cufflinks option (fr-firststrand or fr-secondstrand): ')
            stranded_value = stranded_value.lower().replace(' ','')
            stranded_value = ' --stranded ' + str(stranded_value)
    else:
        stranded_value = ''

    method = input('XPRESSpipe will collect all read quantifications for each sample after processing and save them as a single table.\nWould you like to additionally perform a library normalization on these samples? (yes/no): ')
    if method.lower().replace(' ','') == 'yes':
        method_value = input('Please choose the normalization method (RPM, RPKM, FPKM, TPM): ')
        method_value = method_value.upper().replace(' ','')
        method_value = ' --method ' + str(method_value)
    else:
        method_value = ''

    sjdb = input('Please input the number you provided during reference curation for the --sjdbOverhang option, not what was output in the reference command builder (if no option was used, input \'no\'): ')
    sjdb_value = ' --sjdbOverhang ' + str(sjdb)

    genome = input('Did you use the --genome_size option during reference creation? (yes/no): ')
    if genome.lower().replace(' ','') == 'yes':
        genome_size = input('Please provide that number (integer required): ')
        genome_value = ' --genome_size ' + str(genome_size)
    else:
        genome_value = ''

    processors = input('Would you like to throttle number of processors (default will use all available)? (yes/no): ')
    if processors.lower() == 'no':
        proc_value = ''
    else:
        proc_value = input('What is the maximum number of processors you would like to use? (input must be integer): ')
        proc_value = ' --max_processors ' + str(proc_value)

    # Summarize
    print('\nHere is the command you should use to process your data:\n')
    print(
        'xpresspipe '
        + str(pipeline)
        + ' -i ' + str(indir)
        + ' -o ' + str(outdir)
        + ' -e ' + str(experiment)
        + ' -r ' + str(reference_dir)
        + ' -g ' + str(gtf)
        + str(adaptors)
        + str(twopass)
        + str(quality)
        + str(length)
        + str(rrna_alignments)
        + str(front)
        + str(umi_loc)
        + str(umi_len)
        + str(umi_spacer)
        + str(multimappers)
        + str(duplicates)
        + str(bed)
        + str(quant_value)
        + str(feature_value)
        + str(stranded_value)
        + str(method_value)
        + str(sjdb_value)
        + str(genome_value)
        + str(proc_value)
        )
    print('Please double-check the output for accuracy and that all desired arguments were caught by the builder and incorporated into the command.')

    # Run
    run = input('\nWould you like to run this command now? (yes/no): ')
    if run.lower().replace(' ','') == 'yes':
        os.system(
            'xpresspipe '
            + str(pipeline)
            + ' -i ' + str(indir)
            + ' -o ' + str(outdir)
            + ' -e ' + str(experiment)
            + ' -r ' + str(reference_dir)
            + ' -g ' + str(gtf)
            + str(adaptors)
            + str(twopass)
            + str(quality)
            + str(length)
            + str(rrna_alignments)
            + str(front)
            + str(umi_loc)
            + str(umi_len)
            + str(umi_spacer)
            + str(multimappers)
            + str(duplicates)
            + str(bed)
            + str(quant_value)
            + str(feature_value)
            + str(stranded_value)
            + str(method_value)
            + str(sjdb_value)
            + str(genome_value)
            + str(proc_value)
            )
    else:
        sys.exit(1)


def build_command():
    welcome()

    option = input('What do you want to do:\n1: Curate reference\n2: Process data\nInput choice: ')
    if option == '1':
        ref_guide()
        build_curation()

    elif option == '2':
        run_guide()
        build_pipeline()

    else:
        raise Exception('Invalid input parameter')
