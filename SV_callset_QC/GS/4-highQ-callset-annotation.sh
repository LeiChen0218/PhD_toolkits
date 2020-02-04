#4-highQ-callset-annotation.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
DIR=$ROOT/data/gs/highQ_anno_vcf3
#mkdir -p $DIR
cd $DIR

VCF=$ROOT/data/gs/sample_qc_vcf2/gs.high_conf.joint_samples.vcf.gz
gstools_dir=$ROOT/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh

cd ..
bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo highQ_anno_vcf3/%J-run_annotation.vcf3.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${anno_scr} ${VCF} highQ_anno_vcf3

cd highQ_anno_vcf3
awk '{print $8}' CopyNumberClass.report.dat | sort | uniq -c
#      1 CNCATEGORY
#  19509 DEL
#  14287 DUP
#   9729 MIXED
#   2438 NA

awk '$8 != "NA"' CopyNumberClass.report.dat \
| zjoin -a annotated.output.vcf.gz -b stdin -p "#" -13 -wa \
| bgzip -c > vcf4.rmnonvar.vcf.gz
tabix -p vcf vcf4.rmnonvar.vcf.gz


## variant per sample count

count_scr=$ROOT/scripts/general/CNtools/cn_count_per_sample_20190613.py
bsub -J count_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo count_gs.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python ${count_scr} -v vcf4.rmnonvar.vcf.gz > count_per_sample.gs.vcf4.txt"

samp_key=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/cramid_to_nid.txt
bcftools reheader vcf4.rmnonvar.vcf.gz -s $samp_key -o vcf4.rmnonvar.nid.vcf.gz

###
# SNP LD
LD_scr=$gstools_dir/tagvariants.report_only.sh
SNP_dir=$ROOT/data/general/snp_for_tagvar

#mkdir TagVar_by_chr
cd TagVar_by_chr
for i in {1..22}
do
#mkdir c$i
SNP=$SNP_dir/c$i.snp.reformat_gt.vcf.gz
bsub -J tv${i}_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo c$i/run_TagVar_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${LD_scr} $DIR/vcf4.rmnonvar.nid.vcf.gz $SNP c$i
done

for i in {1..22}
do
awk -v chr="chr"$i '$2==chr' c$i/TagVariants.summary.dat | cut -f 1,6
done > var.maxR2.txt

# wc -l var.maxR2.txt
#   39660 var.maxR2.txt
# awk '$2 >= 0.5 && $2 !="NA"' var.maxR2.txt | wc -l
#   24570


# IRS
#mkdir irs
irs_scr=${gstools_dir}/irs.report_only.sh

bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs/%J-run_irs.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} $DIR/vcf4.rmnonvar.vcf.gz irs/gs.vcf4.irs.report.dat

# awk '$6>1 && $5!="NA" ' *report.dat | wc -l
#    7261
# awk '$6>1 && $5!="NA" && $5>0.5 ' *report.dat | wc -l
#    121

###
# generate bed file for evaluation 
zcat vcf4.rmnonvar.nid.vcf.gz | vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY}' > gs.highQ.bed

##
# generate dosage file for evaluation 
tabix -H vcf4.rmnonvar.nid.vcf.gz | tail -n 1 | cut -f 3,10-6000 > gs.highQ.dosage
zcat vcf4.rmnonvar.nid.vcf.gz | vawk '{print $3,S$*$CNF}' >> gs.highQ.dosage

#mkdir dosage_by_chr
for i in {1..22} X Y 
do
	tabix -H vcf4.rmnonvar.nid.vcf.gz | tail -n 1 | cut -f 3,10-6000 > dosage_by_chr/gs.highQ.c$i.dosage
	tabix vcf4.rmnonvar.nid.vcf.gz chr$i | vawk '{print $3,S$*$CNF}' >> dosage_by_chr/gs.highQ.c$i.dosage
done

cd dosage_by_chr
for i in {1..22} X Y
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo transpose_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash transpose.sh gs.highQ.c$i.dosage gs.highQ.c$i.t.dosage
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
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $corr_scr dosage_by_chr/gs.highQ.c$i.t.dosage  corr_by_chr/gs.corr.c$i.matrix.txt"
done

#mkdir -p matSpDlite/logs
msd_scr=$scr_dir/matSpDlite_calculation.R
for i in {1..22} X Y 
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo matSpDlite/logs/msd_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $msd_scr corr_by_chr/gs.corr.c$i.matrix.txt matSpDlite/matSpDlite.c$i.txt"
done



