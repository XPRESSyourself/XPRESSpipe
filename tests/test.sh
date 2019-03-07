#!/usr/bin/env bash

#############
#TEST INSTALL
#############
xpresspipe --help >> test.out

##########
#TEST TRIM
##########
#Preliminary test
xpresspipe trim --help >> test.out
#Ribosome profiling/SE tests
xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6 >> test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT >> test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q 20 >> test.out
xpresspipe trim -i riboprof_test -o riboprof_out >> test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q -20 >> test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -m 6000 >> test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a None >> test.out
#PE tests
mkdir pe_test_archive
cp pe_test/* pe_test_archive
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None >> test.out
#Remove a file and test if PE catches
rm pe_test/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None >> test.out
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a POLYx >> test.out

#Reset gzipped files
rn -r pe_test
mv pe_test_archive pe_test

#####################
#TEST CREATEREFERENCE
#####################
#Set environment variables
REFERENCE=test_reference
#Preliminary test
xpresspipe createReference --help >> test.out
#Test partial arguments
xpresspipe createReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 >> test.out
#clean up test data
rm -r $REFERENCE/genome
#Test some error-prone tries

#clean up test data
rm -r $REFERENCE/genome
#Create reference for align tests with SE 50bp reads
xpresspipe createReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 --sjdbOverhang 49 >> test.out

#####################
#TEST CREATEREFERENCE
#####################

##############
#TEST TRUNCATE
##############

##############
#TEST MAKEFLAT
##############

###########
#TEST ALIGN
###########
#Set environment variables
TRIMDIR=riboprof_out/trimmed_fastq
#Preliminary test
xpresspipe align --help >> test.out
#Test all arguments
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 --output_bigwig --output_bed --max_processors 10 >> test.out
#Test some error-prone tries

#Final test to use with next steps
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 >> test.out
#clean up test data
rm -r $REFERENCE/genome
#Create reference for align tests with SE 50bp reads
xpresspipe createReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 --sjdbOverhang 49 >> test.out
#Final PE test to use with next steps
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 >> test.out

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

##########################
#RUN FINAL TEST FOR ERRORS
##########################
[[ $(cat test.out | grep -i "error\|exception" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in testing output"; exit 1; }
