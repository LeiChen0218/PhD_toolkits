#!/bin/sh
VAR_ID=$1
WINDOW_SIZE=$2

echo $VAR_ID 

#exp to test the code

# ROOT=
BASE=${ROOT}/3-Candidate_analysis/data/var_plots
CAND_DATA_DIR=${BASE}/${VAR_ID}/data
F_DIR=${BASE}/${VAR_ID}/flanking_${WINDOW_SIZE}

CHR=`cut -f 1 ${F_DIR}/region.*.bed`
REGION=`awk '{print $1":"$2"-"$3}' ${F_DIR}/region.*.bed`


#SNP_VCF=

tabix -h ${SNP_VCF} ${REGION} | bgzip -c > ${F_DIR}/snp.vcf.gz 
tabix -p vcf ${F_DIR}/snp.vcf.gz

LD_SCRIPT=${ROOT}/3-Candidate_analysis/scripts/7.1-sv_ld_with_snps.py 

python ${LD_SCRIPT} \
-t ${CAND_DATA_DIR}/genotype.continous.list \
-v ${F_DIR}/snp.vcf.gz \
-o ${F_DIR}/ld.snp.continous.txt

