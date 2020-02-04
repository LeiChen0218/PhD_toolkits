# overlap_redundancy.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/1-Callset_QC
lumpy_nonbnd=$ROOT/data/lumpy/highQ_anno_vcf3/highQ.nonbnd.bed ### autosome only 
lumpy_bnd=$ROOT/data/lumpy/highQ_anno_vcf3/highQ.bnd.bed
gs_bed=$ROOT/data/gs/highQ_anno_vcf3/gs.highQ.bed
cnvnator_bed=$ROOT/data/cnvnator/highQ_anno_vcf2/cnvnator.highQ.bed

cd $ROOT/data/general
#mkdir genome_density

cd genome_density
CHR_END_REF=/gscmnt/gc2719/halllab/users/lchen/references/b38/hg38.chrom.sizes

POS=1
WINDOW_SIZE=1000000 # 1MB windows
while read CHR END; do
	bash windows_generator.sh ${CHR} ${POS} ${END} ${WINDOW_SIZE}
done < <(head -n 22 ${CHR_END_REF}) | sed 's/[:-]/\t/g' > b38.1mb.windows.bed

bedtools intersect -a b38.1mb.windows.bed -b $gs_bed -wo |  bedtools groupby -g 1,2,3 -c 7 -o count > gs.1mb.density.bed
bedtools intersect -a b38.1mb.windows.bed -b $gs_bed -wo | sort -k1,3V -k8 |  bedtools groupby -g 1,2,3,8 -c 7 -o count > gs.1mb.by_type.density.bed

bedtools intersect -a b38.1mb.windows.bed -b $cnvnator_bed -wo |  bedtools groupby -g 1,2,3 -c 7 -o count > cnvnator.1mb.density.bed
bedtools intersect -a b38.1mb.windows.bed -b $cnvnator_bed -wo | sort -k1,3V -k8 |  bedtools groupby -g 1,2,3,8 -c 7 -o count > cnvnator.1mb.by_type.density.bed

bedtools intersect -a b38.1mb.windows.bed -b $lumpy_bnd -wo |  bedtools groupby -g 1,2,3 -c 7 -o count > lumpy_bnd.1mb.density.bed
#bedtools intersect -a b38.1mb.windows.bed -b $lumpy_bnd -wo | sort -k1,3V -k8 |  bedtools groupby -g 1,2,3,8 -c 7 -o count > lumpy_bnd.1mb.by_type.density.bed

bedtools intersect -a b38.1mb.windows.bed -b $lumpy_nonbnd -wo |  bedtools groupby -g 1,2,3 -c 7 -o count > lumpy_nonbnd.1mb.density.bed
bedtools intersect -a b38.1mb.windows.bed -b $lumpy_nonbnd -wo | sort -k1,3V -k8 |  bedtools groupby -g 1,2,3,8 -c 7 -o count > lumpy_nonbnd.1mb.by_type.density.bed

awk '{OFS="\t";print $1,$2,$3,"BND",$4}' lumpy_bnd.1mb.density.bed |  cat lumpy_nonbnd.1mb.by_type.density.bed - | sort -k1,3V > lumpy.1mb.by_type.density.bed
bedtools groupby -i lumpy.1mb.by_type.density.bed -g 1,2,3 -c 5 -o sum > lumpy.1mb.density.bed




