# aws s3 cp s3://sra-pub-run-odp/sra/SRR1795425/SRR1795425 . --no-sign-request
# aws s3 cp s3://sra-pub-run-odp/sra/SRR1795426/SRR1795426 . --no-sign-request
# ./sratoolkit.3.1.1-ubuntu64/bin/fastq-dump SRR1795425 
# ./sratoolkit.3.1.1-ubuntu64/bin/fastq-dump SRR1795426

cp -r /mnt/c/Users/jorda/Documents/GitHub/XPRESSpipe/tests/periodicity/* ~/periodicity/

REFERENCE=~/periodicity/references
xpresspipe curateReference -o $REFERENCE -f $REFERENCE/fasta -g $REFERENCE/transcripts.gtf --sjdbOverhang 49 --genome_size 723917970

INPUT=~/periodicity/input
OUTPUT=~/periodicity/output
xpresspipe riboseq \
    -i $INPUT \
    -o $OUTPUT \
    -r $REFERENCE \
    --gtf $REFERENCE/transcripts.gtf \
    -e periodicity_file_test \
    -a CTGTAGGCACCATCAAT \
    --sjdbOverhang 49 \
    --quantification_method htseq \
    --cdna_fasta $REFERENCE/fasta/Homo_sapiens.GRCh38.cds.all.fa
