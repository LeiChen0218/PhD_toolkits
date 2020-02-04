2-gscnqual_tuning.sh

gstools_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh
irs_scr=${gstools_dir}/irs.report_only.sh
vcf1=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/gs/anno_vcf0/vcf1.rmnonvar.vcf.gz

cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/gs
# 1.run irs
mkdir irs_vcf1

bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs_vcf1/%J-run_irs.vcf1.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} ${vcf1} irs_vcf1/IRS_vcf1.report.dat

zcat anno_vcf0/vcf1.rmnonvar.vcf.gz \
| vawk '{print $3,I$GSCNCATEGORY, I$GSCNQUAL}' \
> irs_vcf1/id.svtype.qual.txt

sed 's/_/\t/g' irs_vcf1/id.svtype.qual.txt | cut -f 2-4 \
| paste - irs_vcf1/id.svtype.qual.txt > irs_vcf1/vcf1.bed 

# 2.calculate overlap fraction with 1kg+gnomad
overlapsrc_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/Callset_Overlap
anno_src=${overlapsrc_dir}/get_overlapped_list.sh
#anno_src2=${overlapsrc_dir}/get_overlapped_list_merge_small_variants.sh

bed_file_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/bed_other_callsets
SV1kg_bed=${bed_file_dir}/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.bed.gz
gnomad_b38_bed=${bed_file_dir}/gnomad_v2_sv.sites.b38.short.bed
gs5k_bed=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/gs/irs_vcf1/vcf1.bed 

#mkdir overlapFrac_vcf1
cd overlapFrac_vcf1

bash ${anno_src} ${gs5k_bed} ${SV1kg_bed} gs5k.SV1kg.txt
bash ${anno_src} ${gs5k_bed} ${gnomad_b38_bed} gs5k.gnomad_b38.txt

cd ..
mkdir gscnqual_tuning_vcf1
printf "ID\tSVTYPE\tQUAL\tCHR\tSTART\tEND\tIRS_P\tNPROBES\tIN_SV1KG\tIN_GNOMAD\n" > gscnqual_tuning_vcf1/tuning.table

zjoin -a irs_vcf1/id.svtype.qual.txt -b irs_vcf1/IRS_vcf1.report.dat \
| zjoin -a stdin -b overlapFrac_vcf1/gs5k.SV1kg.txt \
| zjoin -a stdin -b overlapFrac_vcf1/gs5k.gnomad_b38.txt \
| cut -f 1-3,5-9,16,18 >>  gscnqual_tuning_vcf1/tuning.table

##
# plot tuning curves in R with all the infomation
# GSCNQUAL cut-off : 2 for all variant type

zcat anno_vcf0/vcf1.rmnonvar.vcf.gz | vawk --header 'I$GSCNQUAL >= 2' | bgzip -c > gscnqual_tuning_vcf1/vcf2.high_conf.vcf.gz 
tabix -p vcf gscnqual_tuning_vcf1/vcf2.high_conf.vcf.gz 




