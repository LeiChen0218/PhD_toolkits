#!/bin/sh
VAR_ID=$1
echo $VAR_ID


#exp to test the code
#VAR_ID="42333"
#TRAIT="height_combined"

#VAR_ID="CNV_chr16_69982813_69983880"
#TRAIT="Pyr_combined"

#ROOT= 

INVCF=${ROOT}/2-Association_test/data/
PHENO=${INVCF}/general/qt_finnseq_20190620.renormalized.wgs.5k_ids.power80.ped

#############
# Automatically make plots for input candidate SV
# 1. Create directory
mkdir -p ${VAR_ID}/{data,plots,logs}


# 2. Extract information specifically for this variant
# # genotype info
if [[ "$VAR_ID" == CNV* ]]; then
VCF=${INVCF}/gs/gs.highQ.mac3.autosome.dechr.vcf.gz
ASSOC_DIR=${INVCF}/gs/epacts_p80qt
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.discrete.list
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.continous.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$CN}' >> ${VAR_ID}/data/genotype.discrete.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$CNF}' >> ${VAR_ID}/data/genotype.continous.list
elif [[ "$VAR_ID" == chr* ]]; then
VCF=${INVCF}/cnvnator/cnvnator.highQ.mac3.autosome.dechr.ref.vcf.gz
ASSOC_DIR=${INVCF}/cnvnator/epacts_p80qt
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.discrete.list
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.continous.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$CN}' >> ${VAR_ID}/data/genotype.discrete.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$CNF}' >> ${VAR_ID}/data/genotype.continous.list
else
VCF=${INVCF}/lumpy/lumpy.highQ.mac3.autosome.dechr.vcf.gz
ASSOC_DIR=${INVCF}/lumpy/epacts_p80qt
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.discrete.list
tabix -H ${VCF} | tail -n 1 | cut -f 3,10-6000 > ${VAR_ID}/data/genotype.continous.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$GT}' >> ${VAR_ID}/data/genotype.discrete.list
zcat ${VCF} | awk -v var=${VAR_ID} '$3==var' | vawk '{print $3,S$*$AB}' >> ${VAR_ID}/data/genotype.continous.list
fi

/bin/bash transpose.sh ${VAR_ID}/data/genotype.discrete.list ${VAR_ID}/data/genotype.discrete.t.list
/bin/bash transpose.sh ${VAR_ID}/data/genotype.continous.list ${VAR_ID}/data/genotype.continous.t.list

printf "SAMPLE\tGT_DIS\tGT_CONT\n" > ${VAR_ID}/data/genotype.t.list
zjoin -a ${VAR_ID}/data/genotype.discrete.t.list -b ${VAR_ID}/data/genotype.continous.t.list |\
cut -f 1,2,4 | tail -n +2 >> ${VAR_ID}/data/genotype.t.list


# # association results
TRAIT_LIST=${INVCF}/general/renormalized.trait.list
printf "TRAIT\tVAR\tCHR\tPOS\tPVALUE\tBETA\n" >  ${VAR_ID}/data/assoc.results.all_traits.sorted.txt
for trait in `cat ${TRAIT_LIST}`
do
	awk -v t=${trait} -v var=${VAR_ID} '{OFS="\t";if( $1==var) print t,$1,$2,$3,$12,$13}' ${ASSOC_DIR}/${trait}/*.${trait}.*.epacts
done | sort -k5,5g >>  ${VAR_ID}/data/assoc.results.all_traits.sorted.txt









