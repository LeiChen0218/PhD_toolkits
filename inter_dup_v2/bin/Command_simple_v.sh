#Commands.sh
DIR='/gscmnt/gc2719/halllab/users/lchen/inter_dup'

mkdir -p $DIR/outputs/no_filter
# Call interspersed duplication patterns from bedpe file
python $DIR/scripts/inter_dup.py \
-b $DIR/inputs/gtex.lumpy.gs.svscore.low_conf.bedpe \
-o $DIR/outputs/no_filter/gtex.lumpy.interdup.50.txt \
-s 50 \
-r 0.8
#real	1m39.023s
#user	1m38.206s
#sys	0m0.412s

# Extract copy numbers from CNVnator
mkdir -p $DIR/no_filter/copynumber/cnvnator_50/log
mkdir -p $DIR/log/cnvnator_50

cat $DIR/outputs/no_filter/gtex.lumpy.interdup.50.txt | cut -f 1 > $DIR/outputs/no_filter/dup_region.50.txt
echo "exit" >> $DIR/outputs/no_filter/dup_region.50.txt

source /gsc/pkg/root/root/bin/thisroot.sh
#GTEX root file for cnvnator
while read SAMPLE BAM 
do
bomb -m 10 -J $SAMPLE.cnvnator -o $DIR/log/cnvnator_50/$SAMPLE.cnvnator.50.log -e $DIR/log/cnvnator_50/$SAMPLE.cnvnator.50.log \
"cat $DIR/outputs/no_filter/dup_region.50.txt | \
/gscmnt/gc2719/halllab/src/speedseq/bin/cnvnator-multi \
-root $BAM \
-genotype 100 > $DIR/outputs/no_filter/copynumber/cnvnator_50/$SAMPLE.50.txt"
done < $DIR/inputs/bam.list

# Test for copy number evidence (svclassfier) 
python $DIR/scripts/cn_evidence.py \
-i $DIR/inputs/gtex.lumpy.gs.svscore.low_conf.vcf \
-v $DIR/outputs/no_filter/gtex.lumpy.interdup.50.txt \
-c $DIR/outputs/no_filter/copynumber/cnvnator_50 \
-o $DIR/outputs/no_filter/copynumber/gtex.lumpy.interdup.cn.50.all.vcf \
-g $DIR/inputs/gender.txt \
-m 'large_sample'




