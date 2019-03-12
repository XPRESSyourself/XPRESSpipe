#!/usr/bin/env bash

#############
#TEST INSTALL
#############
#Preliminary test
xpresspipe --help >> test.out
#Error checking
[[ $(cat test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in INIT testing output"; exit 1; }
rm test.out

##########
#TEST TRIM
##########
mkdir riboprof_out
#Preliminary test
xpresspipe trim --help >> trim_test.out
#Ribosome profiling/SE tests
xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6 >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q 20 >> trim_test.out
xpresspipe trim -i riboprof_test -o riboprof_out >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -q -20 >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a None >> trim_test.out
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a POLYx >> trim_test.out
#Test some errors
[[ $(xpresspipe trim -i riboprof_test/ -o riboprof_out/ -m 6000 | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { exit 0; }
#PE tests
mkdir pe_out
mkdir pe_test_archive
cp pe_test/* pe_test_archive
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None >> trim_test.out
#Remove a file and test if PE catches
rm pe_test/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
[[ $(xpresspipe trim -i pe_test -o pe_out --min_length 50 -a None None | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { exit 0; }
#Reset gzipped files
rm -r pe_test
mv pe_test_archive pe_test
rm -r pe_out
rm -r riboprof_out
mkdir pe_out
mkdir riboprof_out
#Make files for downstream testing
xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT >> trim_test.out
xpresspipe trim -i pe_test -o pe_out --min_length 50 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC >> trim_test.out
#Error checking
[[ $(cat trim_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in TRIM testing output"; exit 1; }
rm trim_test.out

#####################
#TEST MAKEREFERENCE
#####################
#Set environment variables
REFERENCE=test_reference
#Preliminary test
xpresspipe makeReference --help >> makeref_test.out
#Test partial arguments
xpresspipe makeReference -o $REFERENCE/ -f $REFERENCE/ -g $REFERENCE/transcripts.gtf -m 10 >> makeref_test.out
#clean up test data
rm -r $REFERENCE/genome
#Test some error-prone tries

#clean up test data
rm -r $REFERENCE/genome
#Error checking
[[ $(cat makeref_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in MAKEREFERENCE testing output"; exit 1; }
rm makeref_test.out
##############
#TEST TRUNCATE
##############
#Preliminary test
xpresspipe truncate --help >> truncate_test.out
#Use all
xpresspipe truncate -g $REFERENCE/transcripts.gtf -t 15 -c >> truncate_test.out
#clean up test data
rm $REFERENCE/transcripts_*
#Error checking
[[ $(cat truncate_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in TRUNCATE testing output"; exit 1; }
rm truncate_test.out

##############
#TEST MAKEFLAT
##############
#Preliminary test
xpresspipe makeFlat --help >> flatten_test.out
#Run all options
xpresspipe makeFlat -i $REFERENCE >> flatten_test.out
#clean up test data
rm $REFERENCE/transcripts_*
#Error checking
[[ $(cat flatten_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in MAKEFLAT testing output"; exit 1; }
rm flatten_test.out

#####################
#TEST CURATEREFERENCE
#####################
#Preliminary test
xpresspipe curateReference --help >> curate_test.out
#Get references for downstream
SE_REF=se_reference
PE_REF=pe_reference
cp -r $REFERENCE $SE_REF
cp -r $REFERENCE $PE_REF
xpresspipe curateReference -o $SE_REF/ -f $SE_REF/ -g $SE_REF/transcripts.gtf -t 45 -m 10 --sjdbOverhang 49 >> curate_test.out
xpresspipe curateReference -o $PE_REF/ -f $PE_REF/ -g $PE_REF/transcripts.gtf -m 10 >> curate_test.out
#Error checking
[[ $(cat curate_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in MAKEFLAT testing output"; exit 1; }
rm curate_test.out

###########
#TEST ALIGN
###########
#Set environment variables
SE_TRIM=riboprof_out/trimmed_fastq
PE_TRIM=pe_out/trimmed_fastq
#Preliminary test
xpresspipe align --help >> align_test.out
#Test all arguments
xpresspipe align -i $SE_TRIM -o riboprof_out -t SE -r $SE_REF --sjdbOverhang 49 --output_bigwig --output_bed --max_processors 10 >> align_test.out
xpresspipe align -i $SE_TRIM -o riboprof_out -t SE -r $SE_REF --sjdbOverhang 49 >> align_test.out
#Test some error-prone tries

#clean up test data
rm -r riboprof_out/alignments
#Create reference for align tests with SE 50bp reads
xpresspipe align -i $SE_TRIM -o riboprof_out -t SE -r $SE_REF --sjdbOverhang 49 >> align_test.out
#Final PE test to use with next steps
xpresspipe align -i $PE_TRIM -o pe_out -t PE -r $PE_REF >> align_test.out
#Error checking
[[ $(cat align_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in ALIGN testing output"; exit 1; }
rm align_test.out

###########
#TEST COUNT
###########
#Preliminary test
xpresspipe count --help >> count_test.out
#Run some tests
xpresspipe count -i riboprof_out/alignments/ -o riboprof_out/ -r se_reference/ -e se_test >> count_test.out
xpresspipe count -i riboprof_out/alignments/ -o riboprof_out/ -r se_reference/ >> count_test.out
#Run some error prone tests

#Clean up data
rm -r riboprof_out/counts
#Final counting
xpresspipe count -i riboprof_out/alignments/ -o riboprof_out/ -r se_reference/ -e se_test --count_coding --truncate
xpresspipe count -i pe_out/alignments/ -o pe_out/ -r pe_reference/ -e pe_test
#Error checking
[[ $(cat count_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in COUNT testing output"; exit 1; }
rm count_test.out

###############
#TEST NORMALIZE
###############
#Preliminary test
xpresspipe normalizeMatrix --help >> normalize_test.out
#Run some tests
xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.csv --method RPKM -g se_reference/transcripts_coding_truncated.gtf >> normalize_test.out
xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.csv --method RPM
xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.csv --method FPKM -g se_reference/transcripts_coding_truncated.gtf >> normalize_test.out
xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.csv --method LOG
#Run some error prone tests
[[ $(xpresspipe normalizeMatrix -d riboprof_out/counts/se_test_counts_table.csv --method FPKM | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { exit 0; }
#Error checking
[[ $(cat normalize_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in NORMALIZE testing output"; exit 1; }
rm normalize_test.out

############
#TEST DESEQ2
############
#Preliminary test
xpresspipe diffxpress --help >> diffex_test.out
#Run some tests
xpresspipe diffxpress -d test_r/test_dataset.tsv --sample test_r/sample_info.tsv --design replicate+condition
#Run some error prone tests

#Error checking
[[ $(cat diffex_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in DIFFXPRESS testing output"; exit 1; }
rm diffex_test.out

#########################
#TEST BATCH NORMALIZATION
#########################
#Preliminary test
xpresspipe normalizeMatrix --help >> batch_test.out
#Run some tests
xpresspipe normalizeMatrix -d test_r/test_dataset.tsv --batch test_r/batch_info.tsv
#Run some error prone tests

#Error checking
[[ $(cat batch_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in BATCH testing output"; exit 1; }
rm batch_test.out

##############
#TEST METAGENE
##############
xpresspipe metagene --help >> quality_test.out

xpresspipe metagene -i riboprof_out/alignments/ -o riboprof_out -r se_reference -t DEFAULT -e se_test >> quality_test.out

######################
#TEST READDISTRIBUTION
######################
xpresspipe readDistribution --help >> quality_test.out

xpresspipe readDistribution -i riboprof_out/trimmed_fastq -o riboprof_out -e se_test >> quality_test.out

#################
#TEST PERIODICITY
#################
xpresspipe periodicity --help >> quality_test.out

xpresspipe periodicity -i riboprof_out/alignments/ -o riboprof_out -g se_reference/transcripts.gtf -e se_test >> quality_test.out

#Error checking
[[ $(cat batch_test.out | grep -i "error\|exception\|command not found" | wc -l) -eq 0 ]] || { echo "Errors or exceptions were present in BATCH testing output"; exit 1; }
rm batch_test.out

###############
#TEST RRNAPROBE
###############

##################
#TEST CONVERTNAMES
##################

##############
#TEST SERNASEQ
##############

#Test one file

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
rm *STARtmp
rm -r pe_out
rm -r pe_test_archive
rm -r riboprof_out
rm -r *reference/genome
rm -r *reference/genome
rm -r *reference/transcripts_*
rm -r se_reference
rm -r pe_reference
rm test_r/*batched*
rm test_r/*diffx*
