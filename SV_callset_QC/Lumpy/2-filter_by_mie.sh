# 2-filter_by_mie.sh
cd /gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/data/lumpy
vcf=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/0-Callset_generation/data/lumpy/Allsamples.lumpy.one-bad-coord-filtered.vcf.gz
filter_scr=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC/scripts/lumpy/2a-zfilter_mie.sh

bsub -q ccdg \
	 -M 2000000 \
	 -R 'select[mem>2000] rusage[mem=2000]' \
     -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' \
     -oo "logs/$J-filter_mie.log" \
	 bash ${filter_scr} ${vcf} | bgzip -c > Allsamples.lumpy.one-bad-coord-filtered.mie.vcf.gz
tabix -p vcf Allsamples.lumpy.one-bad-coord-filtered.mie.vcf.gz

zcat Allsamples.lumpy.one-bad-coord-filtered.mie.vcf.gz \
| vawk '{print $3,$7,I$MSQ,I$SVTYPE}' > id.filter.msq.svtype
# double check the variant

zcat Allsamples.lumpy.one-bad-coord-filtered.mie.vcf.gz \
| vawk --header '$7=="PASS"' | bgzip -c > lumpy.mie.pass.vcf.gz

tabix -p vcf lumpy.mie.pass.vcf.gz