#!/bin/sh
VAR_ID=$1
WINDOW_SIZE=$2
echo $VAR_ID $WINDOW_SIZE

#exp to test the code
#VAR_ID="44368"
#WINDOW_SIZE=100

# ROOT=
BASE=${ROOT}/3-Candidate_analysis/data/var_plots
CAND_DATA_DIR=${BASE}/${VAR_ID}/data

SAMP_LIST=${CAND_DATA_DIR}/rep100sample.list
F_DIR=${BASE}/${VAR_ID}/flanking_${WINDOW_SIZE}

samp=`head -n 1 ${SAMP_LIST}`
printf "paste %s/cn/%s.cn" ${F_DIR} ${samp} > ${F_DIR}/combine_samples.sh

for samp in `tail -n +2 ${SAMP_LIST}`
do
	printf " <(cut -f 4 %s/cn/%s.cn)" ${F_DIR} ${samp} >> ${F_DIR}/combine_samples.sh
done 

printf " > %s/rep100sample.cn" ${F_DIR} >> ${F_DIR}/combine_samples.sh

chmod 755 ${F_DIR}/combine_samples.sh
/bin/bash ${F_DIR}/combine_samples.sh

