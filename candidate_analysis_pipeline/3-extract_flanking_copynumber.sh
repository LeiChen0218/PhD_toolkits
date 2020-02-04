#!/bin/sh
CAND_CHR=$1
CAND_POS=$2
CAND_END=$3
VAR_DIR=$4
WINDOW_NUMBER=$5

echo $VAR_ID $WINDOW_NUMBER

#exp to test the code
#VAR_ID="chr1:628901-636500"
#WINDOW_SIZE=100

#VAR_ID="CNV_chr16_69982813_69983880"

ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
BASE=${ROOT}/3-Candidate_analysis/data/var_plots
CAND_DATA_DIR=${BASE}/${VAR_DIR}/data

CAND_SIZE="$(($CAND_END-$CAND_POS))"
echo ${CAND_CHR}":"${CAND_POS}"-"${CAND_END}

## get the flanking region (500kb left and right)
CHR_END_REF=/gscmnt/gc2719/halllab/users/lchen/references/b38/b38_chr_end.dechr.txt

# define window size by variant size
WINDOW_SIZE=$(( CAND_SIZE / 10 ))
if [ "$WINDOW_SIZE" -lt 100 ]
	then 
	WINDOW_SIZE=100
fi

REGION_SIZE=$(( WINDOW_SIZE * WINDOW_NUMBER ))

POS=$(( CAND_POS - REGION_SIZE ))
if [ "$POS" -lt 0 ]
	then
	POS=0
fi

END=$(( CAND_END + REGION_SIZE ))
CHR_END=`awk -v chr=$CAND_CHR '$1==chr' ${CHR_END_REF} | cut -f 2`
if [ "$END" -gt "$CHR_END" ]
	then
	END=${CHR_END}
fi

F_DIR=${BASE}/${VAR_DIR}/flanking_${WINDOW_NUMBER}
mkdir -p $F_DIR
printf ${CAND_CHR}"\t"${POS}"\t"${END}"\n" > $F_DIR/region.${REGION_SIZE}.bed 

## extract CNVnator results
/bin/bash windows_generator.sh ${CAND_CHR} ${POS} ${END} ${WINDOW_SIZE} > ${F_DIR}/region.${WINDOW_NUMBER}bp.windows.list
echo "exit" >>  ${F_DIR}/region.${WINDOW_NUMBER}bp.windows.list
sv_dir=${ROOT}/general_info/finn.sample.all.table

mkdir -p ${F_DIR}/cn/logs
#cd ${F_DIR}
## extract CNVnator results

while read x sdir 
do
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
bsub -a 'docker(halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2)' \
-oo ${F_DIR}/cn/logs/$x.cnvnator.log \
-q ccdg -n 5 -M 50000000 \
-R "rusage[gtmp=10, mem=50000] select[mem>50000]" \
-g /lchen/cnvnator \
/bin/bash ${ROOT}/3-Candidate_analysis/scripts/cnvnator_b38.sh $x $sdir ${F_DIR}/region.${WINDOW_NUMBER}bp.windows.list ${F_DIR}
done < <(cut -f 1,5 ${sv_dir} | zjoin -a stdin -b ${CAND_DATA_DIR}/rep100sample.list -wa)















