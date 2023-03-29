#!/bin/bash

BAM=$1
REF=$2
SONIC=$3
BED=$5
OUTFILE=$4

tardis -i $BAM --ref $REF --sonic $SONIC --out ${OUTFILE}$BED
cat "${OUTFILE}.vcf" | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > "${OUTFILE}.tardis.vcf"
# edit the vcf to add the sample name 
sed -i -E "s/(#CHROM.+FORMAT\t)(.+)/\1TD_\2/" "${OUTFILE}.tardis.vcf"