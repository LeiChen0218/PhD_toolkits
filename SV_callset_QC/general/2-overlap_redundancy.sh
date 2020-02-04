# overlap_redundancy.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
lumpy_nonbnd=$ROOT/data/lumpy/highQ_anno_vcf3/highQ.nonbnd.bed ### autosome only 
gs_bed=$ROOT/data/gs/highQ_anno_vcf3/gs.highQ.bed
cnvnator_bed=$ROOT/data/cnvnator/highQ_anno_vcf2/cnvnator.highQ.bed

gs_dosage=$ROOT/data/gs/highQ_anno_vcf3/dosage_by_chr
lumpy_dosage=$ROOT/data/lumpy/highQ_anno_vcf3/dosage_by_chr
cnvnator_dosage=$ROOT/data/cnvnator/highQ_anno_vcf2/dosage_by_chr

cd $ROOT/data/general
# 1. dosage redundancy 
#paste $ROOT/data/gs/highQ_anno_vcf3/gs.samp.list $ROOT/data/lumpy/highQ_anno_vcf3/lumpy.samp.list $ROOT/data/cnvnator/highQ_anno_vcf2/cnvnator.samp.list | awk '$2!=$3'
#mkdir -p dosage_redund/logs
cd dosage_redund

scr_dir=$ROOT/scripts/general/GTRedundancy
corr_scr=$scr_dir/corr_calculation_all_sv.R

for i in {1..22}
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/corr_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $corr_scr $gs_dosage/gs.highQ.c$i.t.dosage $lumpy_dosage/lumpy.highQ.c$i.t.dosage \
$cnvnator_dosage/cnvnator.highQ.c$i.t.dosage all.corr.c$i.matrix.txt"
done

corr_scr_2=$scr_dir/corr_calculation_lumpy_gs.R
for i in X Y
do
bsub -q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/corr_c$i.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash -c \
"Rscript $corr_scr_2 $gs_dosage/gs.highQ.c$i.t.dosage \
$lumpy_dosage/lumpy.highQ.c$i.t.dosage all.corr.c$i.matrix.txt"
done

cd ..
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
"Rscript $msd_scr dosage_redund/all.corr.c$i.matrix.txt matSpDlite/matSpDlite.c$i.txt"
done

##### spacial redundancy
#mkdir spacial_cluster
cd spacial_cluster
bedtools cluster -i $cnvnator_bed -d 10 > cnvnator.highQ.d10clustered.bed
bedtools groupby -i cnvnator.highQ.d10clustered.bed -g 6 -c 4,5 -o count,freqasc > cnvnator.spacial_cluster.var_count.var_type_count.txt

bedtools cluster -i $gs_bed -d 10 > gs.highQ.d10clustered.bed
bedtools groupby -i gs.highQ.d10clustered.bed -g 6 -c 4,5 -o count,freqasc > gs.spacial_cluster.var_count.var_type_count.txt
grep -v "chrX" gs.highQ.d10clustered.bed | grep -v "chrY" > gs.highQ.autosome.d10clustered.bed
bedtools groupby -i gs.highQ.autosome.d10clustered.bed -g 6 -c 4,5 -o count,freqasc > gs.spacial_cluster.autosome.var_count.var_type_count.txt

bedtools cluster -i $lumpy_nonbnd -d 10 > lumpy.highQ.d10clustered.bed
bedtools groupby -i lumpy.highQ.d10clustered.bed -g 6 -c 4,5 -o count,freqasc > lumpy.spacial_cluster.var_count.var_type_count.txt

awk '$5=="DEL" || $5=="DUP"' $lumpy_nonbnd > lumpy.cnv.bed  
bedtools cluster -i lumpy.cnv.bed   -d 10 > lumpy.highQ.cnv.d10clustered.bed
bedtools groupby -i lumpy.highQ.cnv.d10clustered.bed -g 6 -c 4,5 -o count,freqasc > lumpy.cnv.spacial_cluster.var_count.var_type_count.txt

count_scr=$ROOT/scripts/general/CNtools/count_cluster_per_sample.py
gs_vcf=$ROOT/data/gs/highQ_anno_vcf3/vcf4.rmnonvar.nid.vcf.gz
lumpy_vcf=$ROOT/data/lumpy/highQ_anno_vcf3/vcf4.del.dup.gt_only.cram_id.vcf.gz
cnvnator_vcf=$ROOT/data/cnvnator/highQ_anno_vcf2/vcf3.rmnonvar.vcf.gz

bsub -N -u leichen@wustl.edu \
-J ccount_lumpy \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/cluster.per_sample_count.lumpy_cnv.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python $count_scr -v $lumpy_vcf \
-c lumpy.highQ.cnv.d10clustered.bed \
> ClusterSampleCount.lumpy.cnv.txt"

bsub -N -u leichen@wustl.edu \
-J ccount_gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/cluster.per_sample_count.gs.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python $count_scr -v $gs_vcf \
-c gs.highQ.autosome.d10clustered.bed \
-g CN \
> ClusterSampleCount.gs.cnv.txt"

bsub -N -u leichen@wustl.edu \
-J ccount_cnvnator \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo logs/cluster.per_sample_count.cnvnator.log \
-M 1000000 \
-R 'select[mem>1000] rusage[mem=1000]' \
bash -c "python $count_scr -v $cnvnator_vcf \
-c cnvnator.highQ.d10clustered.bed \
-g CN \
> ClusterSampleCount.cnvnator.txt"

###
# callset overlap
cd ..
#mkdir callset_overlap
cd callset_overlap
sv1kg=$ROOT/data/general/bed_other_callsets/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.bed.gz
gnomad=$ROOT/data/general/bed_other_callsets/gnomad_v2_sv.cnvs.b38.short.bed

bedtools intersect -a ../spacial_cluster/lumpy.cnv.bed -b $sv1kg -r -f 0.5 -wo > lumpy.cnv.sv1kg.bed
bedtools intersect -a ../spacial_cluster/lumpy.cnv.bed -b $gnomad -r -f 0.5 -wo > lumpy.cnv.gnomad.bed

bedtools intersect -a $gs_bed -b $sv1kg -r -f 0.5 -wo > gs.cnv.sv1kg.bed
bedtools intersect -a $gs_bed -b $gnomad -r -f 0.5 -wo > gs.cnv.gnomad.bed

bedtools intersect -a $cnvnator_bed -b $sv1kg -r -f 0.5 -wo > cnvnator.cnv.sv1kg.bed
bedtools intersect -a $cnvnator_bed -b $gnomad -r -f 0.5 -wo > cnvnator.cnv.gnomad.bed

bedtools intersect -a ../spacial_cluster/lumpy.cnv.bed -b $gs_bed -r -f 0.5 -wo > lumpy.cnv.gs.bed
bedtools intersect -a ../spacial_cluster/lumpy.cnv.bed -b $cnvnator_bed -r -f 0.5 -wo > lumpy.cnv.cnvnator.bed
bedtools intersect -a $gs_bed -b $cnvnator_bed -r -f 0.5 -wo > gs.cnv.cnvnator.bed





