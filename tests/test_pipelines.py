import os
import sys
import pandas as pd
import scipy.stats as stats

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/scripts/XPRESSyourself/XPRESSpipe/tests/'

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))

# Check outputs produced a file
def check_file(file, type, size=1000000):

    if os.path.isfile(file) == False:
        print('Expected file:', file)
        list_files(__path__)
        print(os.system('tree ' + str(__path__)))
        raise Exception(str(type) + ' file cannot be found')
    else:
        if os.path.getsize(file) < size:
            raise Exception(str(type) + ' file does not look like it produced expected amount of output\nexpected > ' + str(size) + ' bytes, got ' + str(os.path.getsize(file)) + ' bytes.')
        else:
            pass

# Curate references
print('Creating riboseq reference for testing...')
rp_reference = str(__path__) + 'riboseq/rp_reference/'

os.system(
    'xpresspipe curateReference'
    + ' -o ' + str(rp_reference)
    + ' -f ' + str(rp_reference)
    + ' -g ' + str(rp_reference) + 'transcripts.gtf'
    + ' --sjdbOverhang 49'
    + ' -p -t')

# Test riboseq pipeline on test data
print('Running riboseq test pipeline...')
rp_input = str(__path__) + 'riboseq/riboseq_test/'
rp_output = str(__path__) + 'riboseq/riboseq_out/'
if os.path.exists(rp_output):
    os.system(
        'rm -rf ' + str(rp_output))

os.makedirs(rp_output)
rp_gtf = str(__path__) + 'riboseq/rp_reference/transcripts_CT.gtf'
os.system(
    'xpresspipe riboseq' \
    + ' -i ' + str(rp_input) \
    + ' -o ' + str(rp_output) \
    + ' -r ' + str(rp_reference) \
    + ' -g ' + str(rp_gtf) \
    + ' -e riboseq_test' \
    + ' -a CTGTAGGCACCATCAAT' \
    + ' --method RPM' \
    + ' --sjdbOverhang 49')

# Check outputs
check_file(rp_gtf, 'Truncated GTF')

trim1 = rp_output + 'trimmed_fastq/trimmed_SRR2075925_ribo_DMSO_rep1_small.fastq'
trim2 = rp_output + 'trimmed_fastq/trimmed_SRR2075930_rna_DMSO_rep1_small.fastq'
check_file(trim1, 'Trimmed FASTQ')
check_file(trim2, 'Trimmed FASTQ')

align1 = rp_output + 'alignments_coordinates/SRR2075925_ribo_DMSO_rep1_small_dedupMarked.bam'
align2 = rp_output + 'alignments_coordinates/SRR2075930_rna_DMSO_rep1_small_Aligned.sort.bam'
check_file(align1, 'Aligned BAM')
check_file(align2, 'Aligned BAM')

trans_align = rp_output + 'alignments_transcriptome/SRR2075925_ribo_DMSO_rep1_small_Aligned.toTranscriptome.out.bam'
check_file(trans_align, 'Transcriptome Aligned BAM', size = 100000)

count_table = rp_output + 'counts/riboseq_test_count_table_rpmNormalized.tsv'
check_file(count_table, 'Normalized count table', size = 10000)

complexity = rp_output + 'complexity/riboseq_test_library_complexity_1_summary.pdf'
check_file(complexity, 'Complexity summary', size = 10000)

metagene = rp_output + 'metagene/riboseq_test_CDS_metagene_1_summary.pdf'
check_file(metagene, 'Complexity summary', size = 10000)

multiqc = rp_output + 'riboseq_test_multiqc_report.html'
check_file(multiqc, 'Multiqc summary', size = 10000)

if os.path.isdir(rp_output + 'periodicity') == False:
    raise Exception('Periodicity directory cannot be found')

log = rp_output + 'riboseq_test.log'
check_file(log, 'Log', size = 10000)

print('riboseq tests complete')


# Paired-end tests
print('\n\nCreating paired-end reference for testing...')
pe_reference = str(__path__) + 'paired_end/pe_reference/'

os.system(
    'xpresspipe curateReference'
    + ' -o ' + str(pe_reference)
    + ' -f ' + str(pe_reference)
    + ' -g ' + str(pe_reference) + 'transcripts.gtf'
    + ' --genome_size 11'
    + ' -p')

# Test paired end pipeline on test data
print('Running paired-end test pipeline...')
pe_input = str(__path__) + 'paired_end/pe_test/'
pe_output = str(__path__) + 'paired_end/pe_out/'
if os.path.exists(pe_output):
    os.system(
        'rm -rf ' + str(pe_output))

os.makedirs(pe_output)
pe_gtf = str(__path__) + 'paired_end/pe_reference/transcripts_C.gtf'
os.system(
    'xpresspipe peRNAseq' \
    + ' -i ' + str(pe_input) \
    + ' -o ' + str(pe_output) \
    + ' -r ' + str(pe_reference) \
    + ' -g ' + str(pe_gtf) \
    + ' -e pe_test' \
    + ' -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' \
    + ' --method RPM' \
    + ' --sjdbOverhang 100')

# Check outputs
check_file(pe_gtf, 'Protein coding only GTF')

trim1 = pe_output + 'trimmed_fastq/trimmed_UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq'
trim1_1 = pe_output + 'trimmed_fastq/trimmed_UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq'
trim2 = pe_output + 'trimmed_fastq/trimmed_UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq'
trim2_1 = pe_output + 'trimmed_fastq/trimmed_UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq'
check_file(trim1, 'Trimmed FASTQ')
check_file(trim1_1, 'Trimmed FASTQ')
check_file(trim2, 'Trimmed FASTQ')
check_file(trim2_1, 'Trimmed FASTQ')

align1 = pe_output + 'alignments_coordinates/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-ch2._dedupMarked.bam'
align2 = pe_output + 'alignments_coordinates/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-ch2._dedupRemoved.bam'
check_file(align1, 'Aligned BAM')
check_file(align2, 'Aligned BAM')

trans_align = pe_output + 'alignments_transcriptome/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-ch2._Aligned.toTranscriptome.out.bam'
check_file(trans_align, 'Transcriptome Aligned BAM', size = 100000)

count_table = pe_output + 'counts/pe_test_count_table_rpmNormalized.tsv'
check_file(count_table, 'Normalized count table', size = 10000)

complexity = pe_output + 'complexity/pe_test_library_complexity_1_summary.pdf'
check_file(complexity, 'Complexity summary', size = 10000)

metagene = pe_output + 'metagene/pe_test_exon_metagene_1_summary.pdf'
check_file(metagene, 'Complexity summary', size = 10000)

multiqc = pe_output + 'pe_test_multiqc_report.html'
check_file(multiqc, 'Multiqc summary', size = 10000)

if os.path.isdir(pe_output + 'periodicity') == True:
    raise Exception('Periodicity directory should not be created during paired-end processing')

log = pe_output + 'pe_test.log'
check_file(log, 'Log', size = 10000)

print('paired-end tests complete')
