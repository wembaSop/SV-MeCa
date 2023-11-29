#!/bin/bash
set -Eeuo pipefail
BAM=$1
REF=$2
BED=$3
OUTFILE=$4

# read the sample name from bam
SAMPLE=$(samtools samples $BAM|head -n 1|cut -f1)
echo "delly call"$BED" -t ALL -g $REF -o delly_output.bcf $BAM"
delly call$BED -t ALL -g $REF -o delly_output.bcf $BAM
bcftools view -Ov -o $OUTFILE delly_output.bcf 

# edit the vcf to add the sample name 
sed -i -E "s/(#CHROM.+FORMAT\t).+/\1DL_${SAMPLE}/" $OUTFILE