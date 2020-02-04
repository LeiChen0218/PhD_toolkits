# 1-tune_by_mean_sep.sh
rare_vcf=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/0-Callset_generation/data/cnvnator/cnvnator.rare.vcf.gz
common_vcf=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/0-Callset_generation/data/cnvnator/cnvnator.common.vcf.gz

gstools_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/GStoolkit
anno_scr=${gstools_dir}/annotate_svtype_count.vcf_out.sh
irs_scr=${gstools_dir}/irs.report_only.sh

cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/cnvnator

samp_key=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/general_info/nid_to_cramid.txt

bcftools reheader $rare_vcf -s $samp_key -o cnvnator.rare.cramid.vcf.gz
bcftools reheader $common_vcf -s $samp_key -o cnvnator.common.cramid.vcf.gz
tabix -p vcf cnvnator.rare.cramid.vcf.gz
tabix -p vcf cnvnator.common.cramid.vcf.gz

#mkdir irs_rare0 irs_common0

bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs_rare0/%J-run_irs.rare_vcf.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} cnvnator.rare.cramid.vcf.gz irs_rare0/irs_rare0.report.dat

bsub -N -u leichen@wustl.edu \
-J irs-gs \
-q ccdg \
-a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
-oo irs_common0/%J-run_irs.common_vcf.log \
-M 16000000 \
-R 'select[mem>16000] rusage[mem=16000]' \
bash ${irs_scr} cnvnator.common.cramid.vcf.gz irs_common0/irs_common0.report.dat

zcat cnvnator.rare.cramid.vcf.gz \
| vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY,I$DIST,I$NNONREF}' \
> irs_rare.vcf0.bed

zcat cnvnator.common.cramid.vcf.gz \
| vawk '{print $1,$2,I$END,$3,I$GSCNCATEGORY,I$MEANSEP,I$DIPP}' \
> irs_common.vcf0.bed

# 2.calculate overlap fraction with 1kg+gnomad
overlapsrc_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/general/Callset_Overlap
anno_src=${overlapsrc_dir}/get_overlapped_list.sh

bed_file_dir=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/general/bed_other_callsets
SV1kg_bed=${bed_file_dir}/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.bed.gz
gnomad_b38_bed=${bed_file_dir}/gnomad_v2_sv.sites.b38.short.bed
rare_bed=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/cnvnator/irs_rare.vcf0.bed

#mkdir overlapFrac_rare0 
cd overlapFrac_rare0

bash ${anno_src} ${rare_bed} ${SV1kg_bed} rare.SV1kg.txt
bash ${anno_src} ${rare_bed} ${gnomad_b38_bed} rare.gnomad_b38.txt

cd ..
#mkdir rare_tuning
printf "CHR\tSTART\tEND\tID\tSVTYPE\tDIST\tNNONREF\tIRS_P\tNPROBES\tIN_SV1KG\tIN_GNOMAD\n" > rare_tuning/tuning.table

zjoin -a $rare_bed -b irs_rare0/irs_rare0.report.dat -14 | cut -f 1-7,12,13 \
| zjoin -a stdin -b overlapFrac_rare0/rare.SV1kg.txt -14 \
| zjoin -a stdin -b overlapFrac_rare0/rare.gnomad_b38.txt -14 \
| cut -f 1-9,11,13 >>  rare_tuning/tuning.table

common_bed=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/cnvnator/irs_common.vcf0.bed

#mkdir overlapFrac_common0 
cd overlapFrac_common0

bash ${anno_src} ${common_bed} ${SV1kg_bed} common.SV1kg.txt
bash ${anno_src} ${common_bed} ${gnomad_b38_bed} common.gnomad_b38.txt

cd ..
#mkdir common_tuning
printf "CHR\tSTART\tEND\tID\tSVTYPE\tMEAN_SEP\tDIP_P\tIRS_P\tNPROBES\tIN_SV1KG\tIN_GNOMAD\n" > common_tuning/tuning.table

zjoin -a $common_bed -b irs_common0/irs_common0.report.dat -14 | cut -f 1-7,12,13 \
| zjoin -a stdin -b overlapFrac_common0/common.SV1kg.txt -14 \
| zjoin -a stdin -b overlapFrac_common0/common.gnomad_b38.txt -14 \
| cut -f 1-9,11,13 >>  common_tuning/tuning.table

########
# filter sites:
# common: mean_sep >= 0.47 for dup and del, keep all the mixed
# rare: nonref >= 2 for dup, >=5 for del, and >=7 for mixed (aimed at FDR<=0.1)
# rare -- 0.09355976, common -- 0.1402116 

zcat $rare_vcf \
| vawk --header '(I$NNONREF > 1 && I$GSCNCATEGORY=="DUP") || \
                 (I$NNONREF > 4 && I$GSCNCATEGORY=="DEL") || \
                 (I$NNONREF > 6 && I$GSCNCATEGORY=="MIXED")' \
| bgzip -c > rare_tuning/rare.highQ.vcf.gz

zcat $common_vcf \
| vawk --header 'I$GSCNCATEGORY=="MIXED" || I$MEANSEP >= 0.47' \
| bgzip -c > common_tuning/common.highQ.vcf.gz

tabix -p vcf rare_tuning/rare.highQ.vcf.gz
tabix -p vcf common_tuning/common.highQ.vcf.gz

########
# combine the highQ calls
# mkdir combine_common_rare
# check if the sample orders are the same 
# tabix -H common_tuning/common.highQ.vcf.gz | tail -n 1 > temp
# tabix -H rare_tuning/rare.highQ.vcf.gz | tail -n 1 >> temp
# bash transpose.sh temp | awk '$1!=$2'

# deal with ambiguous duplicates
zcat rare_tuning/rare.highQ.vcf.gz | grep -v "##" | cut -f 3,8 \
| zjoin -a common_tuning/common.highQ.vcf.gz -b stdin -p "#" -13 \
| grep -v "#" | cut -f 3,8,4990 > combine_common_rare/duplicates.list



cd combine_common_rare

sed 's/;/\t/g' duplicates.list | cut -f 1,4 | sed 's/DIPP=//g' | awk '$2>0.05' > remove_from_common.txt
sed 's/;/\t/g' duplicates.list | cut -f 1,4 | sed 's/DIPP=//g' | awk '$2<=0.05' > remove_from_rare.txt


tabix -H ../common_tuning/common.highQ.vcf.gz > header 

zcat ../rare_tuning/rare.highQ.vcf.gz | grep -v "#" \
| zjoin -a stdin -b remove_from_rare.txt -v -13 \
| sed 's/END=/MEANSEP=NA;DIPP=NA;NCOMP=NA;END=/g' > cnvnator.highQ.all.vcf

zcat ../common_tuning/common.highQ.vcf.gz | grep -v "#" \
| zjoin -a stdin -b remove_from_common.txt -v -13 \
| sed 's/END=/CN_MED=NA;CN_MAD=NA;DIST=NA;END=/g' >> cnvnator.highQ.all.vcf

sort -k1,1V -k2,2n -k3,3V cnvnator.highQ.all.vcf | cat header - \
| bgzip -c >  cnvnator.highQ.all.sorted.vcf.gz
tabix -p vcf cnvnator.highQ.all.sorted.vcf.gz

# 53396 variants
rm temp
rm header
rm cnvnator.highQ.all.vcf


