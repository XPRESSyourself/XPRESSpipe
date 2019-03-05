#!/usr/bin/env bash

#############
#TEST INSTALL
#############
xpresspipe --help

##########
#TEST TRIM
##########
#Preliminary test
xpresspipe trim --help
#Ribosome profiling/SE tests
xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q 20
xpresspipe trim -i riboprof_test -o riboprof_out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q -20
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -m 6000
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a None
#PE tests
mkdir pe_test_archive
cp pe_test/* pe_test_archive
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None
#Remove a file and test if PE catches
rm pe_test/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a POLYx

#Reset gzipped files
rn -r pe_test
mv pe_test_archive pe_test

###########
#TEST ALIGN
###########

###########
#TEST ALIGN
###########

###########
#TEST COUNT
###########

###############
#TEST NORMALIZE
###############

##############
#TEST METAGENE
##############

######################
#TEST READDISTRIBUTION
######################

#################
#TEST PERIODICITY
#################

##############
#TEST TRUNCATE
##############

##############
#TEST MAKEFLAT
##############
#Combine with make reference?

#####################
#TEST CREATEREFERENCE
#####################

###############
#TEST RRNAPROBE
###############

##################
#TEST CONVERTNAMES
##################

##############
#TEST SERNASEQ
##############

##############
#TEST PERNASEQ
##############

##############
#TEST RIBOPROF
##############
