# 2019-07-12 WES replications for ~5k samples, all three callset

WES_DIR=/gscmnt/gc2802/halllab/lchen/finmetseq_WES/xhmm_outputs_201705_BatchByReadLen  ## use data from this directory to avoid duplicated samples between WES batches
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
DIR=${ROOT}/2-Association_test/data/

LUMPY=$DIR/lumpy/lumpy.highQ.mac3.autosome.dechr.vcf.gz
GS=$DIR/gs/gs.highQ.mac3.autosome.dechr.vcf.gz
CNVNATOR=$DIR/cnvnator/cnvnator.highQ.mac3.autosome.dechr.ref.vcf.gz

###########
# 
cd ${DIR}/lumpy
NONBND=${ROOT}/1-Callset_QC/data/lumpy/highQ_anno_vcf3/highQ.nonbnd.bed
zjoin -a $NONBND -b lumpy.candidate.p3.txt -14 -22 | cut -f 1-6 > lumpy.candidate.p3.nonbnd.bed 

tabix -H $LUMPY | tail -n 1 | cut -f 3,10-6000 > lumpy.candidate.p3.dosage
zjoin -a $LUMPY -b lumpy.candidate.p3.txt -13 -22 -p "#" \
| vawk '{print $3,S$*$AB}' >> lumpy.candidate.p3.dosage

cd ../gs
zjoin -a CopyNumberClass.report.dat -b gs.candidate.p3.txt -22 | cut -f 1,8,10,11 | sed 's/CNV_//1' \
| sed 's/_/\t/1' | sed 's/_/\t/1' | tail -n +2 | awk '{OFS="\t"; print $1,$2,$3,$6,$4,$5}' > gs.candidate.p3.bed 

tabix -H $GS | tail -n 1 | cut -f 3,10-6000 > gs.candidate.p3.dosage
zjoin -a $GS -b gs.candidate.p3.txt -13 -22 -p "#" \
| vawk '{print $3,S$*$CN}' >> gs.candidate.p3.dosage

cd ../cnvnator
zjoin -a CopyNumberClass.report.dat -b cnvnator.candidate.p3.txt -22 | cut -f 1,8,10,11 | sed 's/CNV_//1' \
| sed 's/_/\t/1' | sed 's/_/\t/1' | tail -n +2 | awk '{OFS="\t"; print $1,$2,$3,$6,$4,$5}' > cnvnator.candidate.p3.bed 

tabix -H $CNVNATOR | tail -n 1 | cut -f 3,10-6000 > cnvnator.candidate.p3.dosage
zjoin -a $CNVNATOR -b cnvnator.candidate.p3.txt -13 -22 -p "#" \
| vawk '{print $3,S$*$CNF}' >> cnvnator.candidate.p3.dosage


cd ..
cat */*p3*.bed | sort -k1,1V -k2,2n -k3,3n > general/candidate.p3.bed 
EPACTS_RESULTS=${DIR}/general/candidate.p3.bed

WORK_DIR=${DIR}/general/WES_rep_rm_batch_dups
#mkdir -p ${WORK_DIR}/{inputs,outputs,logs,info}

cd ${WORK_DIR}/inputs
mkdir bams1 bams2

cd ../info
mkdir regions samples targets traits

cd traits
cut -f 6 ${EPACTS_RESULTS} | sed 's/_rn//g' | sort | uniq > ${WORK_DIR}/info/traits/traits.to-run.list 

cd ../regions
cp ${EPACTS_RESULTS} candidate.p_e3.bed

bedtools merge -i candidate.p_e3.bed -c 4,5 -o distinct,distinct

head -n 1  ${DIR}/gs/gs.candidate.p3.dosage > candidate.p3.dosage
tail -n +2  ${DIR}/gs/gs.candidate.p3.dosage >> candidate.p3.dosage
tail -n +2  ${DIR}/cnvnator/cnvnator.candidate.p3.dosage >> candidate.p3.dosage
tail -n +2  ${DIR}/lumpy/lumpy.candidate.p3.dosage >> candidate.p3.dosage
###
# pull in local pc and liftover the coordinates to build37 (for intersecting WES data)

#Successfully converted 4114 records: View Conversions
#Conversion failed on 223 records.    Display failure file    Explain failure messages
#leis-macbook-pro-2:data leichen$ 
#awk '$1!~/chr[1-9A-Za-z]+_/' candidate.p_e3.b37.bed > candidate.p_e3.b37.autosome.bed

awk '{OFS="\t";print $1,$2-5000,$3+5000,$4,$5}' candidate.p_e3.b37.autosome.bed \
| awk '{OFS="\t"; if($2<0) print $1,0,$3,$4,$5; else print $0}' \
| sort -k1,1V -k2,2n -k3,3n > candidate.p_e3.f5kb.b37.bed 

zjoin -a candidate.p_e3.f5kb.b37.bed -b $EPACTS_RESULTS -14 -24 | cut -f 1-5,11 > candidate.p_e3.f5kb.b37.traits.bed

zjoin -a candidate.p_e3.b37.autosome.bed -b $EPACTS_RESULTS -14 -24 | cut -f 1-5,12 | sed 's/^chr//g' > candidate.p_e3.b37.traits.bed

#bedtools merge -i candidate.p_e3.f5kb.b37.bed -c 4,5 -o distinct,distinct \
#| sed 's/^chr//g' > candidate.p_e3.f5kb.b37.merged.bed

bedtools merge -i candidate.p_e3.f5kb.b37.traits.bed -c 4,5,6 -o distinct,distinct,distinct \
| sed 's/^chr//g' > candidate.p_e3.f5kb.b37.merged.traits.bed


# [leichen@blade18-1-7 regions]$ wc -l candidate.p3.f1kb.b37.merged.bed
# 2012 candidate.p3.f1kb.b37.merged.bed

#######
# intersect with WES targets
CNV_COORDINATE=${WORK_DIR}/info/regions/candidate.p_e3.f5kb.b37.merged.traits.bed
BAM_BATCH=1
REGION_BED=${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/regions.bed
bedtools intersect -a ${REGION_BED} -b $CNV_COORDINATE -wo | awk '{OFS="\t"; print $1":"$2"-"$3, $0}' >  ${WORK_DIR}/info/targets/exons.cnvs.intersect.bams${BAM_BATCH}.b37.txt 
BAM_BATCH=2
REGION_BED=${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/regions.bed
bedtools intersect -a ${REGION_BED} -b $CNV_COORDINATE -wo | awk '{OFS="\t"; print $1":"$2"-"$3, $0}' >  ${WORK_DIR}/info/targets/exons.cnvs.intersect.bams${BAM_BATCH}.b37.txt 

####
# extract exons RD data
BAM_BATCH=1
RD_DIR=${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/by_chr_transposed
#RD file name: DATA.PCA_normalized.filtered.sample_zscores.c$CHR.RD.t.txt.gz
# row names: region, col names: samples

# transpose the original WES read depth matrix with all samples
# for i in {1..22}
# do
# bsub -q ccdg \
# -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
# -oo ${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/log/transpose_c$i.log \
# -M 96000000 \
# -R 'select[mem>96000] rusage[mem=96000]' \
# bash transpose.sh ${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/by_chr/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.txt \
# ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.t.txt
# done

# for i in {1..22}
# do
# mv  ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.txt ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.rm_2k_wgs.RD.txt
# bgzip ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.rm_2k_wgs.RD.txt
# #bgzip ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.t.txt
# done

#$ zcat $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c1.RD.t.txt.gz | head -n 1 > ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage
head -n 1 $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c1.RD.t.txt > ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage

for CHR in {1..22}
do
	zjoin -a $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c${CHR}.RD.t.txt \
	      -b ${WORK_DIR}/info/targets/exons.cnvs.intersect.bams${BAM_BATCH}.b37.txt \
	      -wa >> ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage
done 
#  exons: 20147 samples: 9536

# for i in {1..22}
# do
# bgzip ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.t.txt
# done

BAM_BATCH=2
RD_DIR=${WES_DIR}/filtered_zscore/bams${BAM_BATCH}/by_chr_transposed
#RD file name: DATA.PCA_normalized.filtered.sample_zscores.c$CHR.RD.t.txt.gz
# row names: region, col names: samples
#zcat $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c1.RD.t.txt.gz | head -n 1 > ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage

head -n 1 $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c1.RD.t.txt > ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage

for CHR in {1..22}
do
  zjoin -a $RD_DIR/DATA.PCA_normalized.filtered.sample_zscores.c${CHR}.RD.t.txt \
        -b ${WORK_DIR}/info/targets/exons.cnvs.intersect.bams${BAM_BATCH}.b37.txt \
        -wa >> ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.dosage
done 
#  exons: 20139 samples: 9863

# for i in {1..22}
# do
# bgzip ${RD_DIR}/DATA.PCA_normalized.filtered.sample_zscores.c${i}.RD.t.txt
# done


#######
#Make psedo vcf file for epacts to run
CONVERT_SCRIPT=$ROOT/2-Association_test/scripts/general/cn2vcf.py

BAM_BATCH=1
cd ${WORK_DIR}/inputs/bams${BAM_BATCH}/
python $CONVERT_SCRIPT -i finn.5k.cands.rd.${BAM_BATCH}.dosage -o finn.5k.cands.rd.${BAM_BATCH}.vcf
bgzip finn.5k.cands.rd.${BAM_BATCH}.vcf
tabix -p vcf finn.5k.cands.rd.${BAM_BATCH}.vcf.gz

BAM_BATCH=2
cd ${WORK_DIR}/inputs/bams${BAM_BATCH}/
python $CONVERT_SCRIPT -i finn.5k.cands.rd.${BAM_BATCH}.dosage -o finn.5k.cands.rd.${BAM_BATCH}.vcf
bgzip finn.5k.cands.rd.${BAM_BATCH}.vcf
tabix -p vcf finn.5k.cands.rd.${BAM_BATCH}.vcf.gz

#############
# Exclude WGS samples in the test

SAMPLE_TABLE=$ROOT/general_info/finn.samp.all.wcramdir.seqdate.table
zjoin -a ${SAMPLE_TABLE} -b $ROOT/general_info/wgs_samples.list  >  ${WORK_DIR}/info/samples/samples-in-5k-test.all.table
zjoin -a ${SAMPLE_TABLE} -b $ROOT/general_info/wgs_samples.list  | cut -f 2 >  ${WORK_DIR}/info/samples/samples-in-5k-test.oid.list

BAM_BATCH=1
cd ${WORK_DIR}/inputs/bams${BAM_BATCH}/
vcftools --remove ${WORK_DIR}/info/samples/samples-in-5k-test.oid.list \
         --gzvcf finn.5k.cands.rd.${BAM_BATCH}.vcf.gz \
         --out finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp \
         --recode-INFO-all --recode
bgzip finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf
tabix -p vcf finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf.gz
# After filtering, kept 7246 out of 9536 Individuals

BAM_BATCH=2
cd ${WORK_DIR}/inputs/bams${BAM_BATCH}/
vcftools --remove ${WORK_DIR}/info/samples/samples-in-5k-test.oid.list \
         --gzvcf finn.5k.cands.rd.${BAM_BATCH}.vcf.gz \
         --out finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp \
         --recode-INFO-all --recode
bgzip finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf
tabix -p vcf finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf.gz
# After filtering, kept 8166 out of 9863 Individuals

###########################
# Run EMMAX on WES data

epacts=/gscmnt/gc2719/halllab/users/lganel/src/EPACTS-3.2.9/bin/epacts  #epacts for big sample size

pheno_file=/gscmnt/gc2719/halllab/users/lchen/finmetseq_WGS_201610/inputs/qt_finnseq_20161017.ped
kin_mat=/gscmnt/gc2802/halllab/lganel/FinMetSeqExomes/kinship.kinf

export LD_LIBRARY_PATH=/gscmnt/gc2719/halllab/src/gcc-4.9.2/lib64:$LD_LIBRARY_PATH

BAM_BATCH=1
mkdir -p ${WORK_DIR}/outputs/bam${BAM_BATCH}

while read trait
do
mkdir -p ${WORK_DIR}/outputs/bam${BAM_BATCH}/$trait
bsub -n 4 -M 50000000 \
-R 'select[mem>50000] rusage[mem=50000]' \
-g /lchen/SVassoc \
-oo ${WORK_DIR}/logs/wes.emmax.$trait.bam${BAM_BATCH}.log \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-q ccdg \
"$epacts single \
--vcf ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf.gz \
--ped $pheno_file \
--field CN \
--min-maf -1000000 \
--min-mac -1000000 \
--kin $kin_mat \
--pheno $trait \
--test q.emmax \
--out ${WORK_DIR}/outputs/bam${BAM_BATCH}/$trait/finmetseq.wes20k.$trait.emmax.kin.lm.assoc \
--run 5"
done <  ${WORK_DIR}/info/traits/traits.to-run.list

ll ${WORK_DIR}/outputs/bam${BAM_BATCH}/*/*epacts.gz

# cd  ${WORK_DIR}/outputs/bam${BAM_BATCH}

# printf "TRAIT\t" > wes.p3.bam${BAM_BATCH}.epacts
# zcat Alb_combined/*kin.lm.assoc.epacts.gz | head -n 1 >> wes.p3.bam${BAM_BATCH}.epacts
# for d in `ls -d *combined `; 
# do zcat $d/*kin.lm.assoc.epacts.gz |\
# awk -v t=$d '{OFS="\t";if($11 <= 0.001) print t,$0}' | sed 's/\///' \
# >> wes.p3.bam${BAM_BATCH}.epacts; done


BAM_BATCH=2
mkdir -p ${WORK_DIR}/outputs/bam${BAM_BATCH}

while read trait
do
mkdir -p ${WORK_DIR}/outputs/bam${BAM_BATCH}/$trait
bsub -n 4 -M 50000000 \
-R 'select[mem>50000] rusage[mem=50000]' \
-g /lchen/SVassoc \
-oo ${WORK_DIR}/logs/wes.emmax.$trait.bam${BAM_BATCH}.log \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-q ccdg \
"$epacts single \
--vcf ${WORK_DIR}/inputs/bams${BAM_BATCH}/finn.5k.cands.rd.${BAM_BATCH}.rm_wgs_samp.recode.vcf.gz \
--ped $pheno_file \
--field CN \
--min-maf -1000000 \
--min-mac -1000000 \
--kin $kin_mat \
--pheno $trait \
--test q.emmax \
--out ${WORK_DIR}/outputs/bam${BAM_BATCH}/$trait/finmetseq.wes20k.$trait.emmax.kin.lm.assoc \
--run 5"
done <  ${WORK_DIR}/info/traits/traits.to-run.list

# cd  ${WORK_DIR}/outputs/bam${BAM_BATCH}

# printf "TRAIT\t" > wes.p3.bam${BAM_BATCH}.epacts
# zcat Alb_combined/*kin.lm.assoc.epacts.gz | head -n 1 >> wes.p3.bam${BAM_BATCH}.epacts
# for d in `ls -d *combined `; 
# do zcat $d/*kin.lm.assoc.epacts.gz |\
# awk -v t=$d '{OFS="\t";if($11 <= 0.001) print t,$0}' | sed 's/\///' \
# >> wes.p3.bam${BAM_BATCH}.epacts; done


######################
# Meta-analysis of two batches of WES data (BETA,SD) # NOT PVALUE!!!

cd ${WORK_DIR}/outputs
mkdir batch_combined meta

while read trait
do
  zjoin -a bam1/${trait}/*kin.lm.assoc.epacts.gz \
        -b bam2/${trait}/*kin.lm.assoc.epacts.gz \
        -14 -24 |\
        awk '{OFS="\t";print $4,$12,$13,$26,$27}' |\
        tail -n +2 > batch_combined/finn.5k.cands.rd.wes.${trait}.txt
done < ${WORK_DIR}/info/traits/traits.to-run.list


metasoft=/gscmnt/gc2719/halllab/users/lchen/bin/metasoft/Metasoft.jar
ptable=/gscmnt/gc2719/halllab/users/lchen/bin/metasoft/HanEskinPvalueTable.txt
while read trait
do
  /gapp/x64linux/opt/java/jdk/jdk1.8.0_60/bin/java -jar $metasoft \
-input batch_combined/finn.5k.cands.rd.wes.${trait}.txt \
-output meta/finmetseq.${trait}.wes.meta.txt \
-log ${WORK_DIR}/logs/finmetseq.${trait}.meta.log \
-pvalue_table $ptable 
done < ${WORK_DIR}/info/traits/traits.to-run.list

##############
# match back to WGS results
cd ${WORK_DIR}/info/targets

zjoin -wa -a exons.cnvs.intersect.bams1.b37.txt -b exons.cnvs.intersect.bams2.b37.txt | sort | uniq > exons.cnvs.intersect.both.txt
cut -f 1,10 exons.cnvs.intersect.both.txt | sort | uniq > exons.trait.list
cut -f 1,8 exons.cnvs.intersect.both.txt | sort | uniq > exons.cnv.list

while read region traits; do
  T=`echo $traits | sed 's/,/ /g' `
  for t in $T
  do
    echo $region $t
  done
done < exons.trait.list > exons.trait.sep.list

while read region cnvs; do
  T=`echo $cnvs | sed 's/,/ /g' `
  for t in $T
  do
    echo $region $t
  done
done < exons.cnv.list > exons.cnv.sep.list

sed 's/ /\t/g' exons.cnv.sep.list | sort -k2,2V  | bedtools groupby -g 2 -c 1 -o distinct > cnv.exon.table

sed 's/_rn//g'  exons.trait.sep.list > exons.trait.sep.ori.list


####
# choose the target with the highest genotype correlation with the wgs data
scr=$ROOT/2-Association_test/scripts/general/5a-pick_target.py
samp_dic=$ROOT/general_info/finn.samp.all.wcramdir.seqdate.table

bsub -n 4 -M 50000000 \
-R 'select[mem>50000] rusage[mem=50000]' \
-N -u leichen@wustl.edu \
-g /lchen/SVassoc \
-oo logs/cnv.exon.r2.b1.log \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-q ccdg \
bash -c \
"python $scr -l cnv.exon.table \
-v ../regions/candidate.p3.dosage \
-e ${WORK_DIR}/inputs/bams1/finn.5k.cands.rd.1.dosage \
-s $samp_dic -o cnv.exon.bams1.r2.table"

bsub -n 4 -M 50000000 \
-N -u leichen@wustl.edu \
-R 'select[mem>50000] rusage[mem=50000]' \
-g /lchen/SVassoc \
-oo logs/cnv.exon.r2.b2.log \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-q ccdg \
bash -c \
"python $scr -l cnv.exon.table \
-v ../regions/candidate.p3.dosage \
-e ${WORK_DIR}/inputs/bams2/finn.5k.cands.rd.2.dosage \
-s $samp_dic -o cnv.exon.bams2.r2.table"

awk '$3>0.1' cnv.exon.bams1.r2.table | sort -r -k1,1V -k3,3g  \
| bedtools groupby -g 1 -c 2,3 -o first,first > cnv.top_exon.bams1.table

awk '$3>0.1' cnv.exon.bams1.r2.table > cnv.exon.r_01.bams1.table

awk '$3>0.1' cnv.exon.bams2.r2.table | sort -r -k1,1V -k3,3g  \
| bedtools groupby -g 1 -c 2,3 -o first,first > cnv.top_exon.bams2.table

awk '$3>0.1' cnv.exon.bams1.r2.table > cnv.exon.r_01.bams1.table
awk '$3>0.1' cnv.exon.bams2.r2.table > cnv.exon.r_01.bams2.table

# cnv.exon.bams1.r2.table | sort -k3,3g
# CNV_chr6_32036849_32039003  6:32008444-32008549 0.000
# CNV_chr6_32036849_32039003  6:32011525-32011670 0.000
# CNV_chr6_32036849_32039003  6:32008645-32008912 0.001
# CNV_chr6_32036849_32039003  6:32010231-32010384 0.004
# CNV_chr6_32036849_32039003  6:32008182-32008362 0.007
# CNV_chr6_32036849_32039003  6:32009547-32009712 0.007
# CNV_chr6_32036849_32039003  6:32010041-32010141 0.011
# CNV_chr6_32036849_32039003  6:32009125-32009228 0.023 # 5k WES signal
# CNV_chr6_32036849_32039003  6:32007781-32007983 0.025
# CNV_chr6_32036849_32039003  6:32007322-32007425 0.049
# CNV_chr6_32036849_32039003  6:32006199-32006402 0.064
# CNV_chr6_32036849_32039003  6:32006493-32006594 0.066
# CNV_chr6_32036849_32039003  6:32007132-32007235 0.080
# CNV_chr6_32036849_32039003  6:32007519-32007619 0.097
# CNV_chr6_32036849_32039003  6:32006870-32007026 0.268  # 2k WES signal

zjoin -a cnv.top_exon.bams1.table -b ../regions/candidate.p_e3.bed -24 \
| cut -f 2,9 | sed 's/_rn//g' | sort | uniq > top_exon.bams1.trait.list

zjoin -a cnv.top_exon.bams2.table -b ../regions/candidate.p_e3.bed -24 \
| cut -f 2,9 | sed 's/_rn//g' | sort | uniq > top_exon.bams2.trait.list

zjoin -a cnv.exon.r_01.bams1.table -b ../regions/candidate.p_e3.bed -24 \
| cut -f 2,9 | sed 's/_rn//g' | sort | uniq > r01_exon.bams1.trait.list

zjoin -a cnv.exon.r_01.bams2.table -b ../regions/candidate.p_e3.bed -24 \
| cut -f 2,9 | sed 's/_rn//g' | sort | uniq > r01_exon.bams2.trait.list

cat top_exon.bams1.trait.list top_exon.bams2.trait.list | sort | uniq > top_exon.meta.trait.list

cat r01_exon.bams1.trait.list r01_exon.bams2.trait.list | sort | uniq > r01_exon.meta.trait.list

cd ${WORK_DIR}/outputs
printf "TRAIT\tREGION\t" > cands.traits.top_exon.meta.txt
head -n 1 ${WORK_DIR}/outputs/meta/finmetseq.Pyr_combined.wes.meta.txt | cut -f 1,3-8 >> cands.traits.top_exon.meta.txt

while read region trait; do
  grep ${region} ${WORK_DIR}/outputs/meta/finmetseq.${trait}.wes.meta.txt --color=never |\
  awk -v t=${trait} -v r=$region '{OFS="\t"; print t,r,$1,$3,$4,$5,$6,$7,$8}'
done < ${WORK_DIR}/info/targets/top_exon.meta.trait.list >> cands.traits.top_exon.meta.txt

head -n 1 cands.traits.top_exon.meta.txt > cands.traits.top_exon.meta.p3.txt
awk '$4<0.001 || $7<0.001' cands.traits.top_exon.meta.txt >> cands.traits.top_exon.meta.p3.txt

printf "TRAIT\tREGION\t" > cands.traits.r01_exon.meta.txt
head -n 1 ${WORK_DIR}/outputs/meta/finmetseq.Pyr_combined.wes.meta.txt | cut -f 1,3-8 >> cands.traits.r01_exon.meta.txt

while read region trait; do
  grep ${region} ${WORK_DIR}/outputs/meta/finmetseq.${trait}.wes.meta.txt --color=never |\
  awk -v t=${trait} -v r=$region '{OFS="\t"; print t,r,$1,$3,$4,$5,$6,$7,$8}'
done < ${WORK_DIR}/info/targets/r01_exon.meta.trait.list >> cands.traits.r01_exon.meta.txt

head -n 1 cands.traits.r01_exon.meta.txt > cands.traits.r01_exon.meta.p3.txt
awk '$4<0.001 || $7<0.001' cands.traits.r01_exon.meta.txt >> cands.traits.r01_exon.meta.p3.txt


#######################
# per batch check
BAM_BATCH=1
printf "TRAIT\tREGION\t" > cands.traits.bam${BAM_BATCH}.txt
zcat ${WORK_DIR}/outputs/bam${BAM_BATCH}/Pyr_combined/*epacts.gz | head -n 1 | cut -f 4,11,12 >> cands.traits.bam${BAM_BATCH}.txt

while read region trait; do
  zcat ${WORK_DIR}/outputs/bam${BAM_BATCH}/${trait}/*epacts.gz |\
  grep ${region} --color=never |\
  awk -v t=${trait} -v r=$region '{OFS="\t"; print t,r,$4,$11,$12}'
done < ${WORK_DIR}/info/targets/r01_exon.bams1.trait.list >> cands.traits.bam${BAM_BATCH}.txt

head -n 1 cands.traits.bam${BAM_BATCH}.txt > cands.traits.bam${BAM_BATCH}.p3.txt
awk '$4<0.001' cands.traits.bam${BAM_BATCH}.txt >> cands.traits.bam${BAM_BATCH}.p3.txt


BAM_BATCH=2
printf "TRAIT\tREGION\t" > cands.traits.bam${BAM_BATCH}.txt
zcat ${WORK_DIR}/outputs/bam${BAM_BATCH}/Pyr_combined/*epacts.gz | head -n 1 | cut -f 4,11,12 >> cands.traits.bam${BAM_BATCH}.txt

while read region trait; do
  zcat ${WORK_DIR}/outputs/bam${BAM_BATCH}/${trait}/*epacts.gz |\
  grep ${region} --color=never |\
  awk -v t=${trait} -v r=$region '{OFS="\t"; print t,r,$4,$11,$12}'
done < ${WORK_DIR}/info/targets/r01_exon.bams2.trait.list >> cands.traits.bam${BAM_BATCH}.txt

head -n 1 cands.traits.bam${BAM_BATCH}.txt > cands.traits.bam${BAM_BATCH}.p3.txt
awk '$4<0.001' cands.traits.bam${BAM_BATCH}.txt >> cands.traits.bam${BAM_BATCH}.p3.txt

# match wgs results in R script

##################
# Cut-off WES
#cd ${WORK_DIR}/info/targets

### or combined p > 1.6e-6
tail -n +2 wes_wgs_matched.valid.bam1.txt | cut -f 1,2,17 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5 -o min,max,min > wes_wgs_matched.valid.bam1.merged.txt

tail -n +2 wes_wgs_matched.valid.bam2.txt | cut -f 1,2,17 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5 -o min,max,min > wes_wgs_matched.valid.bam2.merged.txt

tail -n +2 wes_wgs_matched.valid.meta.txt | cut -f 1,2,20,22 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5,6 -o min,max,min,min > wes_wgs_matched.valid.meta.merged.txt

# grep VLDL wes_wgs_matched.meta.txt | grep 6:32 | cut -f 1,2,4,5,14,15,20
# XXL_VLDL_P_combined  6:32006870-32007026 0.0175635 -0.0126732  0.000179  -0.6718 4.29768347148065e-05

printf "TRAIT\tCHR\tPOS\tEND\tCNV\tWGS_P\tCOMBINED_P\n" > wes_wgs_matched.valid.bam1.list 
cut -f 6 wes_wgs_matched.valid.bam1.txt | sed 's/\(.*\)_/\1\t/' | cut -f 2 \
| paste wes_wgs_matched.valid.bam1.txt - | cut -f 6-11,17,18 | tail -n +2 | sort -k3,4 -k7,7g \
| bedtools groupby -g 3,4 -c 5,8,1,6,7 -o min,max,first,first,first >> wes_wgs_matched.valid.bam1.list 


printf "TRAIT\tCHR\tPOS\tEND\tCNV\tWGS_P\tCOMBINED_P\n" > wes_wgs_matched.valid.bam2.list 
cut -f 6 wes_wgs_matched.valid.bam2.txt | sed 's/\(.*\)_/\1\t/' | cut -f 2 \
| paste wes_wgs_matched.valid.bam2.txt - | cut -f 6-11,17,18 | tail -n +2 | sort -k3,4 -k7,7g \
| bedtools groupby -g 3,4 -c 5,8,1,6,7 -o min,max,first,first,first >> wes_wgs_matched.valid.bam2.list 

printf "TRAIT\tCHR\tPOS\tEND\tCNV\tWGS_P\tCOMBINED_P\n" > wes_wgs_matched.valid.meta.list
cut -f 10 wes_wgs_matched.valid.meta.txt | sed 's/\(.*\)_/\1\t/' | cut -f 2 \
| paste wes_wgs_matched.valid.meta.txt - | cut -f 10-14,22,25 | tail -n +2 | sort -k2,3 -k6,6g \
| bedtools groupby -g 2,3 -c 4,7,1,5,6 -o min,max,first,first,first >> wes_wgs_matched.valid.meta.list 

### gather sub-threshold p > 1e-5
tail -n +2 wes_wgs_matched.subthre.bam1.txt | cut -f 1,2,17 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5 -o min,max,min > wes_wgs_matched.subthre.bam1.merged.txt

tail -n +2 wes_wgs_matched.subthre.bam2.txt | cut -f 1,2,17 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5 -o min,max,min > wes_wgs_matched.subthre.bam2.merged.txt

tail -n +2 wes_wgs_matched.subthre.meta.txt | cut -f 1,2,20,22 | sed 's/:/\t/' | sed 's/-/\t/' \
| bedtools groupby -g 1,2 -c 3,4,5,6 -o min,max,min,min > wes_wgs_matched.subthre.meta.merged.txt
