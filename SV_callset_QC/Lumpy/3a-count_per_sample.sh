#!/bin/bash

set -xeo pipefail

INVCF=$1
SAMPLE=$2
OUT=$3

for TYPE in DEL DUP INV BND MEI; do
/gscmnt/gc2719/halllab/bin/bcftools view -s $SAMPLE --no-update $INVCF \
| vawk -v TYPE=$TYPE '{ if (I$SVTYPE==TYPE && ! I$SECONDARY) print $10}' \
| cut -d":" -f 1  | sort | uniq -c \
| awk -v SAMPLE=$SAMPLE -v TYPE=$TYPE 'BEGIN {HOMREF=0; HET=0; HOMALT=0} { if ($2=="0/0") HOMREF=$1; else if ($2=="0/1") HET=$1;else if ($2=="1/1") HOMALT=$1 } END {print SAMPLE,TYPE,HOMREF,HET+HOMALT}' OFS="\t"; 
done   > ${OUT}

