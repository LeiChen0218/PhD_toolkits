#!/bin/sh
F_DIR=$1
sv_dir=$2

num_windows=`wc -l ${F_DIR}/region.*.windows.list | awk '{print $1}'`

for s in `cut -f 1 $sv_dir`;
do 
num_rows=`wc -l ${F_DIR}/cn/$s.cn | awk '{print $1}'`
printf "%s\t%s\t%s\n" $s $num_windows $num_rows 
done > ${F_DIR}/cn_output_row_check.txt

awk '$3 != $4' ${F_DIR}/cn_output_row_check.txt | zjoin -a stdin -b $sv_dir -12 -wb > ${F_DIR}/sample.to.re-run.txt