# process the callset for association study
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
DIR=$ROOT/2-Association_test

#mkdir $DIR
cd $DIR
#mkdir data scripts notes
cd data
#mkdir lumpy gs cnvnator general
#cd ../scripts
#mkdir lumpy gs cnvnator general

## Subset callset (autosome only, samples in phenotype data, then allele frequency)
LUMPY=$ROOT/1-Callset_QC/data/lumpy/highQ_anno_vcf3/vcf4.rmnonvar.autosome.vcf.gz 
GS=$ROOT/1-Callset_QC/data/gs/highQ_anno_vcf3/vcf4.rmnonvar.nid.vcf.gz #including sex chromosomes
CNVNATOR=$ROOT/1-Callset_QC/data/cnvnator/highQ_anno_vcf2/cnvnator.highQ.vcf.gz

sample_list=$DIR/data/general/pheno.samp.list

bcftools view $LUMPY -Oz -S $sample_list > lumpy/lumpy.highQ.pheno_samp.vcf.gz
bcftools view $GS -Oz -S $sample_list > gs/gs.highQ.pheno_samp.vcf.gz
bcftools view $CNVNATOR -Oz -S $sample_list > cnvnator/cnvnator.highQ.pheno_samp.vcf.gz
tabix -p vcf lumpy/lumpy.highQ.pheno_samp.vcf.gz
tabix -p vcf gs/gs.highQ.pheno_samp.vcf.gz
tabix -p vcf cnvnator/cnvnator.highQ.pheno_samp.vcf.gz

# redo the allele frequency annotation
svtools afreq lumpy/lumpy.highQ.pheno_samp.vcf.gz | bgzip -c > lumpy/lumpy.highQ.pheno_samp.af.vcf.gz
tabix -p vcf lumpy/lumpy.highQ.pheno_samp.af.vcf.gz

gstools_dir=$ROOT/1-Callset_QC/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh

bsub -N -u leichen@wustl.edu \
-J anno_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo gs/run_annotation.gs.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${anno_scr} gs/gs.highQ.pheno_samp.vcf.gz gs

bsub -N -u leichen@wustl.edu \
-J anno_cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo cnvnator/run_annotation.cnvnator.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${anno_scr} cnvnator/cnvnator.highQ.pheno_samp.vcf.gz cnvnator

# Filter sites based on allele frequency cut-off
# WES: MAC>=3 
# GTEX: MAF>= 0.05 (148*0.05=7.4)
# Preliminary analysis: MSC >= 10 (.2%, MAF~.1%)
#zcat gs/annotated.output.vcf.gz \
#| vawk --header '$1!="chrX" && $1!="chrY" && I$GSNVARIANT >= 3 && I$GSCALLRATE> 0.5' \
#| bgzip -c > gs/gs.highQ.mac3.autosome.vcf.gz

zcat gs/gs.highQ.mac3.autosome.vcf.gz \
| vawk --header 'I$GSNVARIANT >= 10' \
| bgzip -c > gs/gs.highQ.mac10.autosome.vcf.gz

#zcat cnvnator/annotated.output.vcf.gz \
#| vawk --header 'I$GSNVARIANT >= 3' \
#| bgzip -c > cnvnator/cnvnator.highQ.mac3.autosome.vcf.gz

zcat cnvnator/cnvnator.highQ.mac3.autosome.vcf.gz \
| vawk --header 'I$GSNVARIANT >= 10' \
| bgzip -c > cnvnator/cnvnator.highQ.mac10.autosome.vcf.gz

#zcat lumpy/lumpy.highQ.pheno_samp.af.vcf.gz \
#| vawk --header 'I$AF > 0.0007' \
#| bgzip -c > lumpy/lumpy.highQ.mac3.autosome.vcf.gz

zcat lumpy/lumpy.highQ.mac3.autosome.vcf.gz \
| vawk --header 'I$AF > 0.001' \
| bgzip -c > lumpy/lumpy.highQ.mac10.autosome.vcf.gz

tabix -p vcf gs/gs.highQ.mac10.autosome.vcf.gz
tabix -p vcf cnvnator/cnvnator.highQ.mac10.autosome.vcf.gz
tabix -p vcf lumpy/lumpy.highQ.mac10.autosome.vcf.gz

####
# get rid of "chr" prefix for running epacts (Haven't got to re-run this for MAC10 vcfs -- probably not necessary either?)
zcat gs/gs.highQ.mac10.autosome.vcf.gz | sed 's/^chr//g' | bgzip -c >  gs/gs.highQ.mac10.autosome.dechr.vcf.gz
zcat cnvnator/cnvnator.highQ.mac10.autosome.vcf.gz | sed 's/^chr//g' | bgzip -c >  cnvnator/cnvnator.highQ.mac10.autosome.dechr.vcf.gz
zcat lumpy/lumpy.highQ.mac10.autosome.vcf.gz | sed 's/^chr//g' | bgzip -c >  lumpy/lumpy.highQ.mac10.autosome.dechr.vcf.gz
tabix -p vcf gs/gs.highQ.mac10.autosome.dechr.vcf.gz
tabix -p vcf cnvnator/cnvnator.highQ.mac10.autosome.dechr.vcf.gz
tabix -p vcf lumpy/lumpy.highQ.mac10.autosome.dechr.vcf.gz

# add ref column for running epacts
zcat cnvnator/cnvnator.highQ.mac10.autosome.dechr.vcf.gz \
| vawk --header '{$4="N";$6=100;$7="PASS";print}' | bgzip -c > cnvnator/cnvnator.highQ.mac10.autosome.dechr.ref.vcf.gz
tabix -p vcf cnvnator/cnvnator.highQ.mac10.autosome.dechr.ref.vcf.gz

# callset summary 
#zcat lumpy/lumpy.highQ.mac10.autosome.vcf.gz \
#| vawk '{print $1,$2,$3,I$SVTYPE,I$AF}' > lumpy/lumpy.highQ.mac10.autosome.list
awk '$5>0.001' lumpy/lumpy.highQ.mac3.autosome.list > lumpy/lumpy.highQ.mac10.autosome.list

# callrate > 0.5
awk '$7 > 9 && $2>0.5' gs/CopyNumberClass.report.dat | grep -v "chrX" \
| grep -v "chrY" > gs/gs.highQ.mac10.autosome.list

awk '$7>9' cnvnator/CopyNumberClass.report.dat > cnvnator/cnvnator.highQ.mac10.autosome.list


# subset dosage file for redundancy estimation (Re-ran for MAC10)
#mkdir -p gs/dosage_by_chr
do_dir=$ROOT/1-Callset_QC/data/gs/highQ_anno_vcf3/dosage_by_chr

for i in {1..22}
do
	zjoin -a $do_dir/gs.highQ.c$i.dosage -b gs/gs.highQ.mac10.autosome.list \
	-wa -p ID | bash transpose.sh - gs/dosage_by_chr/gs.highQ.mac10.c$i.t.dosage
done

#mkdir -p cnvnator/dosage_by_chr
do_dir=$ROOT/1-Callset_QC/data/cnvnator/highQ_anno_vcf2/dosage_by_chr

#for i in {1..22}
#do
#	zjoin -a $do_dir/cnvnator.highQ.c$i.dosage -b cnvnator/cnvnator.highQ.mac10.autosome.list \
#	-wa -p ID > cnvnator/dosage_by_chr/cnvnator.highQ.mac10.c$i.dosage
#done

for i in {1..22}
do
	zjoin -a cnvnator/dosage_by_chr/cnvnator.highQ.mac3.c$i.dosage -b cnvnator/cnvnator.highQ.mac10.autosome.list \
	-wa -p ID > cnvnator/dosage_by_chr/cnvnator.highQ.mac10.c$i.dosage
done

for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo cnvnator/dosage_by_chr/logs/transpose_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash transpose.sh cnvnator/dosage_by_chr/cnvnator.highQ.mac10.c$i.dosage cnvnator/dosage_by_chr/cnvnator.highQ.mac10.c$i.t.dosage
done

#mkdir -p lumpy/dosage_by_chr
do_dir=$ROOT/1-Callset_QC/data/lumpy/highQ_anno_vcf3/dosage_by_chr

for i in {1..22}
do
	zjoin -a $do_dir/lumpy.highQ.c$i.dosage -b lumpy/lumpy.highQ.mac10.autosome.list \
	-wa -p ID -23 | bash transpose.sh - lumpy/dosage_by_chr/lumpy.highQ.mac10.c$i.t.dosage
done

#mkdir -p general/dosage_redund/logs

scr_dir=$ROOT/1-Callset_QC/scripts/general/GTRedundancy
corr_scr=$scr_dir/corr_calculation_all_sv.R

for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo general/dosage_redund/logs/corr_c$i.mac10.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $corr_scr gs/dosage_by_chr/gs.highQ.mac10.c$i.t.dosage lumpy/dosage_by_chr/lumpy.highQ.mac10.c$i.t.dosage \
cnvnator/dosage_by_chr/cnvnator.highQ.mac10.c$i.t.dosage general/dosage_redund/all.corr.c$i.mac10.matrix.txt"
done

# 
#mkdir -p general/matSpDlite/logs

msd_scr=$scr_dir/matSpDlite_calculation.R
for i in {1..22} X Y
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo general/matSpDlite/logs/msd_c$i.mac10.log \
-J msd_c${i}_mac10 \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $msd_scr general/dosage_redund/all.corr.c$i.mac10.matrix.txt general/matSpDlite/matSpDlite.c$i.mac10.txt"
done

# pwd=/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/general/matSpDlite
cd general/matSpDlite/
for i in {1..22} X Y
do
sed -n '4p' matSpDlite.c$i.mac10.txt
done 

for i in {1..22} X Y
do
tail matSpDlite.c$i.mac10.txt |  sed -n '3p' 
done 

# 26,495.3 tested 