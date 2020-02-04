# 1-re_annotate_rmnonvar.sh
cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/gs

vcf=/gscmnt/gc2758/analysis/genomestrip-batch-finmetseq-b38-2018/data/big-batch/out/Regenotype1-54/regeno.redund.rmdup.all_chr.vcf.gz
gstools_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh

#mkdir anno_vcf0 logs
bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo annotation/%J-run_annotation.vcf0.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${anno_scr} ${vcf} anno_vcf0

cd anno_vcf0
awk '{print $8}' CopyNumberClass.report.dat | sort | uniq -c
#      1 CNCATEGORY
#  33488 DEL
#  55851 DUP
#  21802 MIXED
# 118470 NA

awk '$8 != "NA"' CopyNumberClass.report.dat \
| zjoin -a annotated.output.vcf.gz -b stdin -p "#" -13 -wa \
| bgzip -c > vcf1.rmnonvar.vcf.gz
tabix -p vcf vcf1.rmnonvar.vcf.gz