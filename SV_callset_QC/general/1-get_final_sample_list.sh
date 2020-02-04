# get_final_sample_list.sh
DIR=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/outliers
cd $DIR

zjoin -a lumpy.sample.table -b gs.sample.table -wa \
| zjoin -a stdin -b cnvnator.sample.table -wa \
| zjoin -a stdin -b cn_outliers_to_exclude.txt -wa -v > joint.sample.table
# maintan the same joint sample list since the only two different outliers are from variant calling control

# 16 more samples with abnormally high copy number across the genome (GS)
#zjoin -a cn_outliers_to_exclude.txt -b cnvnator.sample.table | zjoin -a stdin -b gs.sample.table | zjoin -a stdin -b lumpy.sample.table | cut -f 12 | sort | uniq -c
#      3 Dyslipidemia
#      6 FINRISK
#      4 METSIM
#      1 cohort
#      3 corogene

# cut -f 6 joint.sample.table | sort | uniq -c
#    106 Control
#    349 Dyslipidemia
#    148 Eastern.Fin
#   1058 FINRISK
#   2999 METSIM
#    189 corogene
# Total: 4849 


tail -n +2 joint.sample.table | cut -f 1 > joint.sample.nid.list
tail -n +2 joint.sample.table | cut -f 3 > joint.sample.cram_id.list

cd ../../lumpy/sample_qc_vcf2
bcftools view lumpy.mie.rm_outlier.vcf.gz \
-Oz -S $DIR/joint.sample.nid.list > lumpy.mie.joint_samples.vcf.gz
tabix -p vcf lumpy.mie.joint_samples.vcf.gz

cd ../../gs/sample_qc_vcf2
bcftools view gs.high_conf.rm_outlier.vcf.gz \
-Oz -S $DIR/joint.sample.cram_id.list > gs.high_conf.joint_samples.vcf.gz
tabix -p vcf gs.high_conf.joint_samples.vcf.gz

cd ../../cnvnator/sample_qc_vcf1
vcftools --keep $DIR/joint.sample.nid.list \
         --gzvcf cnvnator.highQ.rm_outlier.recode.vcf.gz \
         --out cnvnator.highQ.joint_samples \
         --recode-INFO-all --recode

bgzip cnvnator.highQ.joint_samples.recode.vcf
tabix -p vcf cnvnator.highQ.joint_samples.recode.vcf.gz


