#!/bin/bash
cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/lumpy

match=/gscmnt/gc2802/halllab/lchen/finmetseq_cnw_assoc_b38_s5000/sample_info/sample_namevs_cram_name.csv
trio_old_id=/gscmnt/gc2802/halllab/dlarson/jira/BIO-2602_PanCCDG_MSQ_recalibration/mie_real/sequenced_trios.sorted.fam

mkdir -p mie/logs
cd mie
### trio data
zjoin -a $trio_old_id -b $match -22 -12 -wb | cut -f 1 > kid.id
zjoin -a $trio_old_id -b $match -22 -13 -wb | cut -f 1 > dad.id
zjoin -a $trio_old_id -b $match -22 -14 -wb | cut -f 1 > mom.id

paste <(cut -f 1 $trio_old_id) kid.id dad.id mom.id > sequenced_trios.newheader.sorted.fam

set -eo pipefail
trio_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1a-trios.sh
# 1-extract trio genotype data and SV info
while read FAM KID DAD MOM; do
        vcf=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/0-Callset_generation/data/lumpy/Allsamples.lumpy.one-bad-coord-filtered.vcf.gz
        #echo $vcf
        trio="${KID}_${DAD}_${MOM}"
        #echo "$trio"
        mkdir -p "$trio"
        bsub -g /lchen/ccdg_jobs \
        -q ccdg \
        -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
        -eo "$trio/mie_counts.err" \
        -oo "$trio/mie_counts.out" \
        bash ${trio_scr} $vcf $KID $DAD $MOM "$trio/mie.counts.gz"
done < sequenced_trios.newheader.sorted.fam

while read FAM KID DAD MOM; do
        trio="${KID}_${DAD}_${MOM}" 
        grep "Successfully" $trio/*out | awk -v t=$trio '{print t,$0}'
done < sequenced_trios.newheader.sorted.fam

while read FAM KID DAD MOM; do
        trio="${KID}_${DAD}_${MOM}" 
        zcat $trio/mie.counts.gz | cut -f 2,3,5 | sort | uniq -c |\
        awk '{OFS="\t"; print $1,$2,$3,$4}' |\
        bedtools groupby -g 2,3 -c 1,1,4 -o first,sum,first | \
        awk -v t=$trio '{OFS="\t"; if($5=="MIE") print t,$1,$2,$3/$4; else print t,$1,$2,1-$3/$4}'
done < sequenced_trios.newheader.sorted.fam > mie_rate_per_trio.txt

while read FAM KID DAD MOM; do
        trio="${KID}_${DAD}_${MOM}" 
        zcat $trio/mie.counts.gz | cut -f 2,3,5 | sort | uniq -c 
done < sequenced_trios.newheader.sorted.fam 


printf "perl /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1b-combine_counts.pl " > /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1c-combine.sh
while read FAM KID DAD MOM; do
        trio="${KID}_${DAD}_${MOM}"
        printf " <(zcat %s/mie.counts.gz)" $trio >> /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1c-combine.sh
done < sequenced_trios.newheader.sorted.fam

bash /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/1c-combine.sh > mie_counts.all_trios.txt
