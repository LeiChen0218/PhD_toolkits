#!/bin/sh
CAND_CHR=$1
CAND_POS=$2
CAND_END=$3
WINDOW_SIZE=$4
WINDOW_NUMBER=$5
OUT_DIR=$6

## get the flanking region 
CHR_END_REF=/gscmnt/gc2719/halllab/users/lchen/references/b38/b38_chr_end.dechr.txt

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

F_DIR=${OUT_DIR}/flanking_${WINDOW_SIZE}
mkdir -p $F_DIR
printf ${CAND_CHR}"\t"${POS}"\t"${END}"\n" > $F_DIR/region.${REGION_SIZE}.bed 


/bin/bash /gscuser/leichen/bin/windows_generator.sh ${CAND_CHR} ${POS} ${END} ${WINDOW_SIZE} > ${F_DIR}/region.${WINDOW_SIZE}bp.windows.list
echo "exit" >>  ${F_DIR}/region.${WINDOW_SIZE}bp.windows.list

