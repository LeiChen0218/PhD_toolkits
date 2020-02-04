#!/bin/sh
dir_list=$1
F_DIR=$2
WINDOW_SIZE=$3

## extract CNVnator results

mkdir -p ${F_DIR}/cn/logs
## extract CNVnator results
scr_cnvnator=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/3-Candidate_analysis/scripts/cnvnator_b38.sh
# use cnvnator_b38.full_name.sh for 1000 GP samples where the full file names are provided

while read x sdir 
do
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
bsub -a 'docker(halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2)' \
-oo ${F_DIR}/cn/logs/$x.cnvnator.log \
-q ccdg -n 5 -M 50000000 \
-R "rusage[gtmp=10, mem=50000] select[mem>50000]" \
-g /lchen/cnvnator \
/bin/bash ${scr_cnvnator} $x $sdir/ ${F_DIR}/region.${WINDOW_SIZE}bp.windows.list ${F_DIR}
done < $dir_list

