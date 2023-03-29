#!/bin/bash

BAM=$1
REF=$2
OUTFILE=$3
MEMORY=$4
CPU=$5
BED=$6

echo "
BAM=$1
REF=$2
OUTFILE=$3
MEMORY=$4
CPU=$5
BED=$6
"
# read the sample name from bam
SAMPLE=$(samtools samples $BAM|head -n 1|cut -f1)

MEM=${MEMORY/ GB/}
echo "configManta.py --bam $BAM --referenceFasta $REF --runDir ./$BED"
configManta.py --bam $BAM --referenceFasta $REF --runDir ./$BED
./runWorkflow.py -j $CPU -g $MEM
cp results/variants/diploidSV.vcf.gz .
bgzip -d diploidSV.vcf.gz
mv diploidSV.vcf $OUTFILE

# edit the vcf to add the sample name 
sed -i -E "s/(#CHROM.+FORMAT\t).+/\1MT_${SAMPLE}/" $OUTFILE