# 3-sample_level_qc.vcf2.sh

INVCF=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/lumpy/lumpy.mie.pass.vcf.gz

cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/lumpy
#mkdir sample_qc_vcf2
cd sample_qc_vcf2

tabix -H $INVCF | tail -n 1 | cut -f 10-6000 | tr '\t' '\n' > vcf2.sample.list

count_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/3a-count_per_sample.sh

#mkdir per_sample

for s in `cat vcf2.sample.list`
do
bsub -J count_lumpy \
-g /lchen/ccdg_jobs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo per_sample/count_$s.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash ${count_scr} ${INVCF} $s per_sample/$s.count 
done

grep "Successfully" per_sample/*log | wc -l
# 5065
cat per_sample/*count > var_count.per_sample.vcf2.txt

#########
# autosome only counts 
zcat $INVCF | awk '$1 != "chrX" && $1 != "chrY"' | bgzip -c > lump.mie.pass.autosome.vcf.gz
tabix -p vcf lump.mie.pass.autosome.vcf.gz

for s in `cat vcf2.sample.list`
do
bsub -J count_lumpy \
-g /lchen/ccdg_jobs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo per_sample/count_$s.autosome.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash ${count_scr} lump.mie.pass.autosome.vcf.gz $s per_sample/$s.autsome.count 
done

grep "Successfully" per_sample/*.autosome.log | wc -l
# 5065
cat per_sample/*.autsome.count > var_count.per_sample.vcf2.autsome.txt

##############
# after plot the samples per variant count, confirmed the artefact in first sequencing batch 
# trying to get rid of those FP calls
# *For LUMPY, the FP sites are all BNDs, so make a specail bed file for BNDs 
zcat lump.mie.pass.autosome.vcf.gz | vawk '{if(I$SVTYPE=="BND") print $1,$2-500,$2+500,$3,"BND"}' > lump.mie.pass.bnd.bed 
zcat lump.mie.pass.autosome.vcf.gz | vawk '{if(I$SVTYPE !="BND") print $1,$2,I$END,$3,I$SVTYPE}' > lump.mie.pass.nonbnd.bed 

cohort15=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.seqdate2015.table
repeats=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/simpleRepeat.b38.sorted.bed

qc_gen_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/
#mkdir -p ${qc_gen_dir}/u1repeats
grep "u=1;" $repeats --color=never | sed 's/^/chr/g' > $qc_gen_dir/u1repeats/u1.repeats.bed 
grep "u=2;" $repeats --color=never | sed 's/^/chr/g' > $qc_gen_dir/u1repeats/u2.repeats.bed # just in case

repeat_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/Callset_Overlap/get_overlapped_list.any_overlap.sh
af_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/CNtools/af_cal_in_subgroup.py
#mkdir rm_cohort15_artefact
cd rm_cohort15_artefact
cut -f 1 $cohort15 | tail -n +2 > nid.2015.list
zcat ../lump.mie.pass.autosome.vcf.gz | python ${af_scr} -s nid.2015.list -g GT > lumpy.mie.pass.autosome.Freq_carrier.2015.txt
cp ../../../cnvnator/sample_qc_vcf1/check_cohort15_artefact/nid.201*.list .

bsub -J count_lumpy \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/lumpy.2016.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"zcat ../lump.mie.pass.autosome.vcf.gz | python ${af_scr} -s nid.2016.list -g GT > lumpy.mie.pass.autosome.Freq_carrier.2016.txt"

bsub -J count_lumpy \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/lumpy.2017.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c \
"zcat ../lump.mie.pass.autosome.vcf.gz | python ${af_scr} -s nid.2017.list -g GT > lumpy.mie.pass.autosome.Freq_carrier.2017.txt"


bash ${repeat_scr} ../lump.mie.pass.bnd.bed $qc_gen_dir/u1repeats/u1.repeats.bed lump.mie.pass.bnd.u1repeat.txt
bash ${repeat_scr} ../lump.mie.pass.bnd.bed $qc_gen_dir/u1repeats/u2.repeats.bed lump.mie.pass.bnd.u2repeat.txt

bash ${repeat_scr} ../lump.mie.pass.nonbnd.bed $qc_gen_dir/u1repeats/u1.repeats.bed lump.mie.pass.nonbnd.u1repeat.txt
bash ${repeat_scr} ../lump.mie.pass.nonbnd.bed $qc_gen_dir/u1repeats/u2.repeats.bed lump.mie.pass.nonbnd.u2repeat.txt

zjoin -b lump.mie.pass.bnd.u1repeat.txt -a lumpy.mie.pass.autosome.Freq_carrier.2015.txt -p "SV_id" \
| zjoin -a stdin -b lump.mie.pass.bnd.u2repeat.txt -p "SV_id"  \
| sed 's/SV_carrier_rate/SV_carrier_rate\tid\tin_u1repeat\tid_2\tin_u2repeat/' > bnd.repeat.Freq_carrier.2015.txt

# just in case 
zjoin -b lump.mie.pass.nonbnd.u1repeat.txt -a lumpy.mie.pass.autosome.Freq_carrier.2015.txt -p "SV_id" \
| zjoin -a stdin -b lump.mie.pass.nonbnd.u2repeat.txt -p "SV_id"  \
| sed 's/SV_carrier_rate/SV_carrier_rate\tid\tin_u1repeat\tid_2\tin_u2repeat/' > nonbnd.repeat.Freq_carrier.2015.txt


######
# outlier defined by fisher's exact test 
cat c*.outlier.NfisherP_200.txt | cut -f 1 | sort | uniq > outlier.list 
# 612 sites
zjoin -a $INVCF -b outlier.list  -v -p "#" -13 | bgzip -c > lumpy.mie.batch_fp_filtered.vcf.gz
tabix -p vcf lumpy.mie.batch_fp_filtered.vcf.gz


####
# another per sample variant count 
zcat lumpy.mie.batch_fp_filtered.vcf.gz | awk '$1 != "chrX" && $1 != "chrY"' | bgzip -c > lumpy.mie.batch_fp_filtered.autosome.vcf.gz
tabix -p vcf lumpy.mie.batch_fp_filtered.autosome.vcf.gz

cd ..

for s in `cat vcf2.sample.list`
do
bsub -J count_lumpy \
-g /lchen/ccdg_jobs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo per_sample/batch_fp_filtered.count_$s.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash ${count_scr} rm_cohort15_artefact/lumpy.mie.batch_fp_filtered.autosome.vcf.gz $s per_sample/batch_fp_filtered.autosome.$s.count 
done

grep "Successfully" per_sample/batch_fp_filtered.*.log | wc -l
# 5065
cat per_sample/batch_fp_filtered.*.count > var_count.per_sample.vcf2.batch_fp_filtered.txt


grep "Successfully" per_sample/cohort15_filtered.*.autosome.log | wc -l
# 5065
cat per_sample/batch_fp_filtered.*.autosome.count > var_count.per_sample.vcf2.cohort15_filtered.autosome.txt

####
# filter lumpy vcf 
cut -f 1 lumpy.var_count.outliers.finns.table | tail -n +2 > samp.to.exclude.list
bcftools view rm_cohort15_artefact/lumpy.mie.batch_fp_filtered.vcf.gz \
-Oz -S ^samp.to.exclude.list > lumpy.mie.rm_outlier.vcf.gz

tabix -p vcf lumpy.mie.rm_outlier.vcf.gz

# combine all the outlier samples and exclude them from all the vcf files
 ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
 cd ${ROOT}/data/general
 #mkdir outliers
lumpy=${ROOT}/data/lumpy/sample_qc_vcf2/lumpy.mie.rm_outlier.vcf.gz
sample_table=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/finn.samp.all.wcramdir.seqdate.table

tabix -H $lumpy | tail -n 1 | cut -f 10-6000 | tr '\t' '\n' \
| zjoin -a $sample_table -b stdin -wa -p "nid" > outliers/lumpy.sample.table #5062





