#3-highQ-callset-annotation.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
DIR=$ROOT/data/lumpy/highQ_anno_vcf3
#mkdir -p $DIR
cd $DIR

VCF=$ROOT/data/lumpy/sample_qc_vcf2/lumpy.mie.joint_samples.vcf.gz

svtools afreq $VCF | bgzip -c > lumpy.highQ.afreq_reanno.vcf.gz
tabix -p vcf lumpy.highQ.afreq_reanno.vcf.gz

zcat lumpy.highQ.afreq_reanno.vcf.gz | vawk --header 'I$AF != 0' | bgzip -c > vcf4.rmnonvar.vcf.gz
tabix -p vcf vcf4.rmnonvar.vcf.gz 

# var per sample count 

tabix -H lumpy.highQ.afreq_reanno.vcf.gz | tail -n 1 | cut -f 10-6000 | tr '\t' '\n' > vcf3.sample.list

####
# exclude small variants
zjoin -a vcf4.rmnonvar.vcf.gz -b ../var_size_filter/small_than_50bp.list -p "#" -v -13 \
| bgzip -c > lumpy.final.highQ.vcf.gz
tabix -p vcf lumpy.final.highQ.vcf.gz

tabix -h lumpy.final.highQ.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | bgzip -c > vcf4.rmnonvar.autosome.vcf.gz
tabix -p vcf vcf4.rmnonvar.autosome.vcf.gz

#mkdir per_sample
count_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/3a-count_per_sample.sh

for s in `cat vcf3.sample.list`
do
bsub -J count_lumpy \
-g /lchen/ccdg_jobs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo per_sample/batch_fp_filtered.count_$s.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash ${count_scr} vcf4.rmnonvar.autosome.vcf.gz $s per_sample/vcf4.autosome.$s.count 
done


grep "Successfully" per_sample/*log | wc -l
# 4865

for s in `cat vcf3.sample.list`
do
	cat per_sample/*$s*count >> var_count.per_sample.autosome.vcf4.txt	
done
####
# run SNP LD 
gstools_dir=$ROOT/scripts/general/GStoolkit
LD_scr=$gstools_dir/tagvariants.report_only.sh

SNP_dir=$ROOT/data/general/snp_for_tagvar

#mkdir TagVar_by_chr
cd TagVar_by_chr
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo simplify_vcf.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"bcftools annotate -x INFO,^FORMAT/GT $DIR/lumpy.final.highQ.vcf.gz \
	| bgzip -c > vcf4.rmnonvar.short.gt_only.vcf.gz"


for i in {1..22}
do
SNP=$SNP_dir/c$i.snp.reformat_gt.vcf.gz
bsub -J tv${i}_lumpy \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo c$i/run_TagVar_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${LD_scr} vcf4.rmnonvar.short.gt_only.vcf.gz $SNP c$i
done

for i in {1..22}
do
awk -v chr="chr"$i '$2==chr' c$i/TagVariants.summary.dat | cut -f 1,6
done > var.maxR2.txt

zjoin -a var.maxR2.txt -b ../../var_size_filter/small_than_50bp.list -v >  var.maxR2.size_filtered.txt

sed 's/_1//g' var.maxR2.size_filtered.txt | sed 's/_2//g' | grep -v "NA" | bedtools groupby -g 1 -c 2 -o max > var.maxR2.bnd_collapsed.txt

# sed 's/_1//g' var.maxR2.size_filtered.txt | sed 's/_2//g' | cut -f 1 | sort | uniq | wc -l
# 35742
# awk '$2>=0.5' var.maxR2.bnd_collapsed.txt | wc -l
# 22645  

# IRS
cd ..
samp_key=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/nid_to_cramid.txt
bcftools reheader $DIR/lumpy.final.highQ.vcf.gz -s $samp_key -o lumpy.vcf4.cramid.vcf.gz
tabix -p vcf lumpy.vcf4.cramid.vcf.gz

BCFTOOL=/gscmnt/gc2719/halllab/bin/bcftools
$BCFTOOL annotate -x ^FORMAT/GT lumpy.vcf4.cramid.vcf.gz \
| vawk --header 'I$SVTYPE == "DEL" || I$SVTYPE == "DUP" ' | bgzip -c \
> vcf4.del.dup.gt_only.cram_id.vcf.gz
tabix -p vcf vcf4.del.dup.gt_only.cram_id.vcf.gz

#mkdir irs
irs_scr=${gstools_dir}/irs.report_only.sh

bsub -N -u leichen@wustl.edu \
-J irs-lumpy \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs/%J-run_irs.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} vcf4.del.dup.gt_only.cram_id.vcf.gz irs/lumpy.vcf4.irs.report.dat

# awk '$6>1 && $5!="NA" ' lumpy.vcf4.irs.report.dat | wc -l
# awk '$6>1 && $5>0.5 && $5!="NA" ' lumpy.vcf4.irs.report.dat | wc -l
# 13*2/3264 = 0.007965686

###
# generate bed file for evaluation 
bnd_pre=$ROOT/data/lumpy/sample_qc_vcf2/lump.mie.pass.bnd.bed
nonbnd_pre=$ROOT/data/lumpy/sample_qc_vcf2/lump.mie.pass.nonbnd.bed

zcat lumpy.final.highQ.vcf.gz | grep -v "#" | cut -f 3 > highQ.var.list
zjoin -a $bnd_pre -b highQ.var.list -14 -wa > highQ.bnd.bed 
zjoin -a $nonbnd_pre -b highQ.var.list -14 -wa > highQ.nonbnd.bed 

# generate dosage file for evaluation 
tabix -H lumpy.final.highQ.vcf.gz | tail -n 1 | cut -f 3,10-6000 > lumpy.highQ.dosage
zcat lumpy.final.highQ.vcf.gz | vawk '{if(! I$SECONDARY) print $3,S$*$AB}' >> lumpy.highQ.dosage

#mkdir dosage_by_chr
for i in {1..22} X Y 
do
	tabix -H lumpy.final.highQ.vcf.gz | tail -n 1 | cut -f 3,10-6000 > dosage_by_chr/lumpy.highQ.c$i.dosage
	tabix lumpy.final.highQ.vcf.gz chr$i | vawk '{if(! I$SECONDARY) print $3,S$*$AB}' >> dosage_by_chr/lumpy.highQ.c$i.dosage
done


cd dosage_by_chr
for i in {1..22} X Y
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo transpose_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash transpose.sh lumpy.highQ.c$i.dosage lumpy.highQ.c$i.t.dosage
done

cd ..
#mkdir -p corr_by_chr/logs

scr_dir=$ROOT/scripts/general/GTRedundancy
corr_scr=$scr_dir/corr_calculation_one_callset.R

for i in {1..22} X Y
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo corr_by_chr/logs/corr_c$i.log \
-M 6000000 \
-R 'select[mem>6000] rusage[mem=6000]' \
bash -c \
"Rscript $corr_scr dosage_by_chr/lumpy.highQ.c$i.t.dosage  corr_by_chr/lumpy.corr.c$i.matrix.txt"
done

#mkdir -p matSpDlite/logs
msd_scr=$scr_dir/matSpDlite_calculation.R
for i in {1..22} X Y 
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo matSpDlite/logs/msd_c$i.log \
-M 16000000 \
-R 'select[mem>6000] rusage[mem=6000]' \
bash -c \
"Rscript $msd_scr corr_by_chr/lumpy.corr.c$i.matrix.txt matSpDlite/matSpDlite.c$i.txt"
done
