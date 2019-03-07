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
xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6 >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q 20 >> trim_test.out
xpresspipe trim -i riboprof_test -o riboprof_out >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q -20 >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -m 6000 >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a None >> trim_test.out
#PE tests
mkdir pe_test_archive
cp pe_test/* pe_test_archive
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None >> trim_test.out
#Remove a file and test if PE catches
rm pe_test/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None >> trim_test.out
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a POLYx >> trim_test.out

#Reset gzipped files
rn -r pe_test
mv pe_test_archive pe_test

[[ $(cat trim_test.out | grep -i "error\|exception" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in TRIM testing output"; exit 1; }
#####################
#TEST MAKEREFERENCE
#####################
#Set environment variables
REFERENCE=test_reference
#Preliminary test
xpresspipe makeReference --help >> makeref_test.out
#Test partial arguments
xpresspipe makeReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 >> makeref_test.out
#clean up test data
rm -r $REFERENCE/genome
#Test some error-prone tries

#clean up test data
rm -r $REFERENCE/genome
#Create reference for align tests with SE 50bp reads
xpresspipe makeReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 --sjdbOverhang 49 >> makeref_test.out

rm -r $REFERENCE/genome

[[ $(cat makeref_test.out | grep -i "error\|exception" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in MAKEREFERENCE testing output"; exit 1; }
##############
#TEST TRUNCATE
##############



##############
#TEST MAKEFLAT
##############


#####################
#TEST CURATEREFERENCE
#####################
SE_REF=se_reference
PE_REF=pe_reference


###########
#TEST ALIGN
###########
#Set environment variables
TRIMDIR=riboprof_out/trimmed_fastq
#Preliminary test
xpresspipe align --help >> align_test.out
#Test all arguments
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 --output_bigwig --output_bed --max_processors 10 >> align_test.out
#Test some error-prone tries

#Final test to use with next steps
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 >> align_test.out
#clean up test data
rm -r $REFERENCE/genome
#Create reference for align tests with SE 50bp reads
xpresspipe createReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -t 10 --sjdbOverhang 49 >> align_test.out
#Final PE test to use with next steps
xpresspipe align -i $TRIMDIR -o riboprof_out -t SE -r $REFERENCE --sjdbOverhang 49 >> align_test.out



[[ $(cat align_test.out | grep -i "error\|exception" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in ALIGN testing output"; exit 1; }
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

###########################
#FINAL CLEANUP OF TEST DATA
###########################
rm .DS_Store
rm *out
rm -r pe_out
rm -r pe_test_archive
rm -r riboprof_out
rm -r *reference/genome
rm -r *reference/genome
