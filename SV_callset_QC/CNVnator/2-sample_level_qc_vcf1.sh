#2-sample_level_qc_vcf1.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
INVCF=${ROOT}/data/cnvnator/combine_common_rare/cnvnator.highQ.all.sorted.vcf.gz
cd $ROOT/data/cnvnator

# mkdir sample_qc_vcf1
cd sample_qc_vcf1

count_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/CNtools/cn_count_per_sample_20190613.py

bsub -J count_cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_cnvnator.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v ${INVCF} > count_per_sample.cnvnator.vcf1.txt"
 
##############
# after plot the samples per variant count, haven't observed the artefact in first sequencing batch 
# still, a sanity check 
zcat $INVCF | vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY}' | bgzip -c > cnvnator.highQ.bed.gz 

#mkdir check_cohort15_artefact
cd check_cohort15_artefact

af_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/CNtools/af_cal_in_subgroup.py
zcat $INVCF | python ${af_scr} -s ${ROOT}/data/lumpy/sample_qc_vcf2/rm_cohort15_artefact/nid.2015.list -g CN > cnvnator.highQ.Freq_carrier.2015.txt
cohort_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general
# observe the enrichment/absence of variants in 2015 cohort, check if it's still due to the simple repeats

qc_gen_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/
cohort15=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.seqdate2015.table
repeat_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/Callset_Overlap/get_overlapped_list.any_overlap.sh

bash ${repeat_scr} ../cnvnator.highQ.bed.gz  $qc_gen_dir/u1repeats/u1.repeats.bed cnvnator.high_conf.u1repeat.txt
bash ${repeat_scr} ../cnvnator.highQ.bed.gz  $qc_gen_dir/u1repeats/u2.repeats.bed cnvnator.high_conf.u2repeat.txt

zjoin -a cnvnator.highQ.Freq_carrier.2015.txt -b cnvnator.high_conf.u1repeat.txt -p "SV_id" \
| zjoin -a stdin -b cnvnator.high_conf.u2repeat.txt -p "SV_id" \
| sed 's/SV_carrier_rate/SV_carrier_rate\tid\tin_u1repeat\tid_2\tin_u2repeat/' > cnvnator.repeat.Freq_carrier.2015.txt

# check other cohort
cut -f 1 $cohort_dir/finn.seqdate2016.table | tail -n +2 > nid.2016.list
cut -f 1 $cohort_dir/finn.seqdate2017.table | tail -n +2 > nid.2017.list
cp ../../../lumpy/sample_qc_vcf2/rm_cohort15_artefact/nid.2015.list .
zcat $INVCF | python ${af_scr} -s nid.2016.list -g CN > cnvnator.highQ.Freq_carrier.2016.txt
zcat $INVCF | python ${af_scr} -s nid.2017.list -g CN > cnvnator.highQ.Freq_carrier.2017.txt

####### observed the artefacts in each cohort, try to do t-test for each subset and get rid of those
# test_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/cnvnator/2c-cnf_ttest_subgroup.py

# for year in 2015 2016 2017
# do
# bsub -q ccdg \
# -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
# -oo logs/t_test.$year.log \
# bash -c \
# "zcat $INVCF | python ${test_scr} -s nid.$year.list -g CNF > cnvnator.highQ.t_test.$year.txt"
# done 

########
# Finnally decided to use Fisher's exact test for enrichment/depletion , and identified ~1k sites for filtering
zjoin -a $INVCF -b outlier_by_NfisherP_200.txt -p '#' -v -13 | bgzip -c > cnvnator.highQ.batch_fp_filtered.vcf.gz
tabix -p vcf cnvnator.highQ.batch_fp_filtered.vcf.gz
#####
# variant per sample count 
cd ..
bsub -J count_cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/count_cnvnator.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.vcf.gz > count_per_sample.cnvnator.batch_fp_filtered.txt"

###
# correct the header:
#tabix -H check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.vcf.gz | tail -n 1 >> vcf_header.txt
#zcat check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.vcf.gz | grep -v "#" \
#| cat vcf_header.txt - | bgzip -c > check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.reheader.vcf.gz
#tabix -p vcf check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.reheader.vcf.gz

cut -f 1 sample_outliers.all_rm.5mad.txt > outlier.list
####
# filter out the samples
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo rm_outlier.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"vcftools --remove outlier.list \
         --gzvcf check_cohort15_artefact/cnvnator.highQ.batch_fp_filtered.vcf.gz \
         --out cnvnator.highQ.rm_outlier \
         --recode-INFO-all --recode"
bgzip cnvnator.highQ.rm_outlier.recode.vcf
tabix -p vcf cnvnator.highQ.rm_outlier.recode.vcf.gz


# combine all the outlier samples and exclude them from all the vcf files
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
 cd ${ROOT}/data/general
 #mkdir outliers
 cnvnator=${ROOT}/data/cnvnator/sample_qc_vcf1/cnvnator.highQ.rm_outlier.recode.vcf.gz
 sample_table=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.samp.all.wcramdir.seqdate.table

tabix -H $cnvnator | tail -n 1 | cut -f 10-6000 | tr '\t' '\n' \
| zjoin -a $sample_table -b stdin -wa -p "nid" > outliers/cnvnator.sample.table 




