#3-sample_level_qc_vcf2.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
INVCF=${ROOT}/data/gs/gscnqual_tuning_vcf1/vcf2.high_conf.vcf.gz 
cd $ROOT/data/gs

#mkdir sample_qc_vcf2
cd sample_qc_vcf2

count_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/CNtools/cn_count_per_sample_20190613.py

#tabix -h $INVCF chr20 | bgzip -c > chr20.text.vcf.gz
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_gs.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v ${INVCF} > count_per_sample.gs.vcf2.txt"
 
samp_dic=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.samp.all.wcramdir.seqdate.table

##############
# after plot the samples per variant count, confirmed the artefact in first sequencing batch 
# trying to get rid of those FP calls

zcat $INVCF | vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY}' > gs.high_conf.bed 

qc_gen_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/
cohort15=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.seqdate2015.table
repeat_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/Callset_Overlap/get_overlapped_list.any_overlap.sh
af_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/CNtools/af_cal_in_subgroup.py

#mkdir rm_cohort15_artefact
cd rm_cohort15_artefact
bash ${repeat_scr} ../gs.high_conf.bed  $qc_gen_dir/u1repeats/u1.repeats.bed gs.high_conf.u1repeat.txt
bash ${repeat_scr} ../gs.high_conf.bed  $qc_gen_dir/u1repeats/u2.repeats.bed gs.high_conf.u2repeat.txt

cut -f 3 $cohort15 | tail -n +2 > cram_id.2015.list
cohort_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general
cut -f 3 $cohort_dir/finn.seqdate2016.table | tail -n +2 > cram_id.2016.list
cut -f 3 $cohort_dir/finn.seqdate2017.table | tail -n +2 > cram_id.2017.list

zcat $INVCF | python ${af_scr} -s cram_id.2015.list -g CN > gs.high_conf.Freq_carrier.2015.txt
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_gs.cohort15_filtered.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"zcat $INVCF | python ${af_scr} -s cram_id.2016.list -g CN > gs.high_conf.Freq_carrier.2016.txt"

bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/cohort17_af.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"zcat $INVCF | python ${af_scr} -s cram_id.2017.list -g CN > gs.high_conf.Freq_carrier.2017.txt"

zjoin -a gs.high_conf.Freq_carrier.2015.txt -b gs.high_conf.u1repeat.txt -p "SV_id" \
| zjoin -a stdin -b gs.high_conf.u2repeat.txt -p "SV_id" \
| sed 's/SV_carrier_rate/SV_carrier_rate\tid\tin_u1repeat\tid_2\tin_u2repeat/' > gs.repeat.Freq_carrier.2015.txt

#######
# site filtered by fisher's exact test
cat c*.outlier.NfisherP_200.txt | cut -f 1 | sort | uniq > outlier.list 
# filter out 740 variants
zjoin -a $INVCF -b outlier.list  -v -p "#" -13 | bgzip -c > gs.high_conf.batch_fp_filtered.vcf.gz
tabix -p vcf gs.high_conf.batch_fp_filtered.vcf.gz

# Similarly to LUMPY, diff(carrier_freq_2015, carrier_freq_non2015) > 0.5 seemed to be a good cut-off
#awk '$2-$3 > 0.5' gs.high_conf.Freq_carrier.2015.txt | cut -f 1 > gs.cohort15.fp.list
#zjoin -a $INVCF -b gs.cohort15.fp.list -v -p "#" -13 | bgzip -c > gs.high_conf.cohort15.vcf.gz
#tabix -p vcf gs.high_conf.cohort15.vcf.gz

cd ..
####
# another per sample variant count 
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_gs.batch_fp_filtered.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v rm_cohort15_artefact/gs.high_conf.batch_fp_filtered.vcf.gz > count_per_sample.gs.batch_fp_filtered.txt"


# check the overlap of outlier samples with previous experiments
#pre_outlier=/gscmnt/gc2802/halllab/lchen/finmetseq_sv_assoc_b38_s5000/info/sample_info/outliers.s72.cramid.list
#outlier=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/gs/sample_qc_vcf2/gs.var_count.outliers.finns.table

#zjoin -a $pre_outlier -b $outlier -v

cut -f 1 gs.var_count.outliers.finns.table | tail -n +2  > samp.to.exclude.list
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo rm_outlier.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"bcftools view -S ^samp.to.exclude.list rm_cohort15_artefact/gs.high_conf.batch_fp_filtered.vcf.gz \
-Oz > gs.high_conf.rm_outlier.vcf.gz"

tabix -p vcf gs.high_conf.rm_outlier.vcf.gz

# combine all the outlier samples and exclude them from all the vcf files
 ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
 cd ${ROOT}/data/general
 #mkdir outliers
 gs=${ROOT}/data/gs/sample_qc_vcf2/gs.high_conf.rm_outlier.vcf.gz
# cnvnator=${ROOT}/data/cnvnator/sample_qc_vcf1/outlier.cram_id.list
 sample_table=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.samp.all.wcramdir.seqdate.table

tabix -H $gs | tail -n 1 | cut -f 10-6000 | tr '\t' '\n' \
| zjoin -a $sample_table -b stdin -wa -13 -p "nid" > outliers/gs.sample.table #4966

 