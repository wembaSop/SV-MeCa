#!/bin/bash

BAM=$1
REF=$2
OUTFILE=$3
PREFIX=$4
CPU=$5
BED="$6"



# read the sample name from bam
SAMPLE=$(samtools samples $BAM|head -n 1|cut -f1)

echo "smoove call -x --genotype --name $PREFIX --outdir . -f $REF --processes ${CPU} $BED $BAM"
smoove call -x --genotype --name $PREFIX --outdir . -f $REF --processes ${CPU} $BED $BAM
bgzip -d *genotyped.vcf.gz
mv *genotyped.vcf $OUTFILE

# edit the vcf to add the sample name 
sed -i -E "s/(#CHROM.+FORMAT\t).+/\1LP_${SAMPLE}/" $OUTFILE