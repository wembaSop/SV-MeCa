#!/bin/bash

BAM=$1
REF=$2
BED=$3
PREFIX=$4

# read the sample name from bam
SAMPLE=$(samtools samples $BAM|head -n 1|cut -f1)

# create the config file for breakdancer max
bam2cfg.pl -g $BAM > "${PREFIX}.cfg"

# run brakdancer max with "-h" to print alle frequency column
breakdancer-max -h "${PREFIX}.cfg" > "${PREFIX}.ctx"

# convert the ctx output to vcf - Only taking DEL and INS in consideration
breakdancertovcf.py -o "${PREFIX}.raw.vcf" $REF "${PREFIX}.ctx"

# Save the header of the vcf
grep "^#" "${PREFIX}.raw.vcf" > "${PREFIX}.unsorted.vcf"

# then add the ID: just enumerate from 1...
grep -v "^#" "${PREFIX}.raw.vcf"| awk 'BEGIN{OFS="\t"};{$3="BD_"NR; print $0}'>> "${PREFIX}.unsorted.vcf"

# sort the vcf file 
cat "${PREFIX}.unsorted.vcf" | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > "${PREFIX}.sorted.vcf"

# compress, index and filter to save regions in the bed file
bgzip "${PREFIX}.sorted.vcf"
tabix -p vcf "${PREFIX}.sorted.vcf.gz"
bcftools view -R $BED -Ov -o "${PREFIX}.vcf" "${PREFIX}.sorted.vcf.gz"

# edit the vcf to add the sample name 
sed -i -E "s/(#CHROM.+FORMAT\t).+/\1BD_${SAMPLE}/" "${PREFIX}.vcf"