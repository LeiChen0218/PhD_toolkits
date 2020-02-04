#3-highQ-callset-annotation.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
DIR=$ROOT/data/cnvnator/highQ_anno_vcf2
#mkdir -p $DIR
cd $DIR

VCF=$ROOT/data/cnvnator/sample_qc_vcf1/cnvnator.highQ.joint_samples.recode.vcf.gz
gstools_dir=$ROOT/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh

cd ..
bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo highQ_anno_vcf2/%J-run_annotation.vcf2.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${anno_scr} ${VCF} highQ_anno_vcf2

cd highQ_anno_vcf2
awk '{print $8}' CopyNumberClass.report.dat | sort | uniq -c
#      1 CNCATEGORY
#  15424 DEL
#  13312 DUP
#  25057 MIXED
#    459 NA

awk '$8 != "NA"' CopyNumberClass.report.dat \
| zjoin -a annotated.output.vcf.gz -b stdin -p "#" -13 -wa \
| bgzip -c > vcf3.rmnonvar.vcf.gz
tabix -p vcf vcf3.rmnonvar.vcf.gz

## variant per sample count

count_scr=$ROOT/scripts/general/CNtools/cn_count_per_sample_20190613.py
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_cnvnator.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v vcf3.rmnonvar.vcf.gz > count_per_sample.cnvnator.vcf3.txt"

### run snp ld 
LD_scr=$gstools_dir/tagvariants.report_only.sh

SNP_dir=$ROOT/data/general/snp_for_tagvar

#mkdir TagVar_by_chr
cd TagVar_by_chr
for i in {1..22}
do
#mkdir c$i
SNP=$SNP_dir/c$i.snp.reformat_gt.vcf.gz
bsub -J tv${i}_cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo c$i/run_TagVar_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${LD_scr} $DIR/vcf3.rmnonvar.vcf.gz $SNP c$i
done

for i in {1..22}
do
awk -v chr="chr"$i '$2==chr' c$i/TagVariants.summary.dat | cut -f 1,6
done > var.maxR2.txt

#  wc -l var.maxR2.txt
#   52090 var.maxR2.txt
# awk '$2 >= 0.5 && $2 !="NA"' var.maxR2.txt | wc -l
#   23991

# IRS
cd ..
samp_key=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/nid_to_cramid.txt
bcftools reheader $DIR/vcf3.rmnonvar.vcf.gz -s $samp_key -o cnvnator.vcf3.cramid.vcf.gz
tabix -p vcf cnvnator.vcf3.cramid.vcf.gz

#mkdir irs
irs_scr=${gstools_dir}/irs.report_only.sh

bsub -N -u leichen@wustl.edu \
-J irs-cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs/%J-run_irs.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} cnvnator.vcf3.cramid.vcf.gz irs/cnvnator.vcf3.irs.report.dat

#awk '$6>1 && $5!="NA" ' *report.dat | wc -l
#    5109
#awk '$6>1 && $5!="NA" && $5>0.5 ' *report.dat | wc -l
#     230

###
# generate bed file for evaluation 
zcat vcf3.rmnonvar.vcf.gz | vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY}' > cnvnator.highQ.bed

##
# generate dosage file for evaluation 
tabix -H vcf3.rmnonvar.vcf.gz | tail -n 1 | cut -f 3,10-6000 > cnvnator.highQ.dosage
zcat vcf3.rmnonvar.vcf.gz | vawk '{print $3,S$*$CNF}' >> cnvnator.highQ.dosage

#mkdir dosage_by_chr
for i in {1..22}
do
	tabix -H vcf3.rmnonvar.vcf.gz | tail -n 1 | cut -f 3,10-6000 > dosage_by_chr/cnvnator.highQ.c$i.dosage
	tabix vcf3.rmnonvar.vcf.gz chr$i | vawk '{print $3,S$*$CNF}' >> dosage_by_chr/cnvnator.highQ.c$i.dosage
done

cd dosage_by_chr
for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo transpose_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash transpose.sh cnvnator.highQ.c$i.dosage cnvnator.highQ.c$i.t.dosage
done

cd ..
#mkdir -p corr_by_chr/logs

scr_dir=$ROOT/scripts/general/GTRedundancy
corr_scr=$scr_dir/corr_calculation_one_callset.R

for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo corr_by_chr/logs/corr_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $corr_scr dosage_by_chr/cnvnator.highQ.c$i.t.dosage  corr_by_chr/cnvnator.corr.c$i.matrix.txt"
done

#mkdir -p matSpDlite/logs
msd_scr=$scr_dir/matSpDlite_calculation.R
for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo matSpDlite/logs/msd_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $msd_scr corr_by_chr/cnvnator.corr.c$i.matrix.txt matSpDlite/matSpDlite.c$i.txt"
done

cp vcf3.rmnonvar.vcf.gz cnvnator.highQ.vcf.gz
tabix -p vcf cnvnator.highQ.vcf.gz



