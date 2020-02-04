#!/bin/bash

set -xeo pipefail

INVCF=$1
# THESE ARE IN FAM/PED file order
KID=$2
DAD=$3
MOM=$4
countfile=$5

# TODO
# Update to calculate MIEs on the fly
#| /gscmnt/gc2719/halllab/bin/bcftools query -f '%TYPE\t%FILTER\t%INFO/RECLASSIFIED\t%INFO/MSQ[\t%GT]\n' \


/gscmnt/gc2719/halllab/bin/bcftools view -s ${KID},${DAD},${MOM} --no-update $INVCF \
| grep -v '^chrY\|^chrX' \
| /gscmnt/gc2719/halllab/bin/bcftools view -g ^miss --no-update - \
| /gscmnt/gc2719/halllab/bin/bcftools query -f '%INFO/SVTYPE\t%FILTER\t%INFO/MSQ[\t%GT]\n' \
| grep "/1"  \
| python /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1ai-classify_mie.py \
| sort -k1,4 \
| bedtools groupby -g 1,2,3,4 -c 1 -o count \
| sed "s/^/${KID}:${DAD}:${MOM}\t/" \
| bgzip -c > $countfile
 



