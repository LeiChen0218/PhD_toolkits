#!/bin/sh
VAR=$1
WINDOW_SIZE=$2
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
BASE=${ROOT}/3-Candidate_analysis/data/var_plots
sv_dir=${ROOT}/general_info/finn.sample.all.table

F_DIR=${BASE}/${VAR}/flanking_${WINDOW_SIZE}

num_windows=`wc -l ${F_DIR}/region.*.windows.list | awk '{print $1}'`

for s in `cat ${BASE}/${VAR}/data/rep100sample.list`;
do 
num_rows=`wc -l ${F_DIR}/cn/$s.cn | awk '{print $1}'`
printf "%s\t%s\t%s\t%s\n" $VAR $s $num_windows $num_rows 
done > ${F_DIR}/cn_output_row_check.txt

awk '$3 != $4' ${F_DIR}/cn_output_row_check.txt | zjoin -a stdin -b $sv_dir -12 -wb > ${F_DIR}/sample.to.re-run.txt