#!/bin/sh
VCF_DIR=$1
TRAIT=$2
echo $VCF $TRAIT

VCF=${VCF_DIR}/snp*vcf.gz

#PHENO_FILE=
#KIN_MAT=
export LD_LIBRARY_PATH=/gscmnt/gc2719/halllab/src/gcc-4.9.2/lib64:$LD_LIBRARY_PATH
EPACTS=/gscmnt/gc2719/halllab/users/lganel/src/EPACTS-3.2.9/bin/epacts 
OUT_DIR=${VCF_DIR}/snp_assoc/${TRAIT}

mkdir -p $OUT_DIR

${EPACTS} single \
--test q.emmax \
--kinf ${KIN_MAT} \
--field PL \
--unit 500000 \
-vcf ${VCF} \
--ped ${PHENO_FILE} \
--out ${OUT_DIR}/snp_assoc_${TRAIT} \
--pheno ${TRAIT} \
--run 1
