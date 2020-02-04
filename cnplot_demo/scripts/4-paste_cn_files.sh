#!/bin/sh
F_DIR=$1
SAMP_LIST=$2

samp=`head -n 1 ${SAMP_LIST} | cut -f 1`
printf "paste %s/cn/%s.cn" ${F_DIR} ${samp} > ${F_DIR}/combine_samples.sh

for samp in `tail -n +2 ${SAMP_LIST} | cut -f 1`
do
	printf " <(cut -f 4 %s/cn/%s.cn)" ${F_DIR} ${samp} >> ${F_DIR}/combine_samples.sh
done 

printf " > %s/all_sample.cn" ${F_DIR} >> ${F_DIR}/combine_samples.sh

chmod 755 ${F_DIR}/combine_samples.sh
/bin/bash ${F_DIR}/combine_samples.sh

