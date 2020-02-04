# interspersed duplication example
# pwd=/gscmnt/gc2719/halllab/users/lchen/projects/inter_dup/example

interdup=/gscmnt/gc2719/halllab/users/lchen/projects/inter_dup/scripts/inter_dup.py
vcf=gtex.lumpy.gs.svscore.low_conf.vcf

# ----
# 1. convert vcf file to bedpe by svtools
# ----
svtools vcftobedpe -i $vcf -o gtex.lumpy.gs.svscore.low_conf.bedpe

# ----
# 2. run inter_dup.py to identify interspersed duplication by pattern search
# ----
python $interdup -b gtex.lumpy.gs.svscore.low_conf.bedpe \
-o gtex.lumpy.gs.svscore.low_conf.txt \
-s 50 \
-r 0.8 

#-b gtex.lumpy.gs.svscore.low_conf.bedpe  #input bedpe file
#-o gtex.lumpy.gs.svscore.low_conf.txt  #output candidates file 
#-s 50   #allow 50bp distance between two "overlapped" breakend, default=0
#-r 0.8  #minium pearson r to define two BNDs as the same event, default=0

# Output file format (GT from the BND with higher dosage was kept as GT for the candidate):
#Column names	Descriptions	Examples
#DUP_pos	Duplication position	1:2811951-2812177
#INS_pos	Insertion position	1:2811132
#BND_1	Source BND	LUMPY_BND_184411
#BND_2	Source BND	LUMPY_BND_184412
#TYPE	Insertion orientation	Inverted
#stat	Pearson's R of two BNDs	0.845687501656
#GT_table	Genotype agreement table (numer of samples : without two BNDs; only with BND1; only with BND2; with both BNDs)	145;1;0;1
#FORMAT	Similar to #format in vcf	GT:AB
#GTEX-N7MS (SAMPLE column )	GT:AB value for each sample (max in BNDs)	0/0:0
#…		
#(Other SAMPLE columns)		
#…		


# ----
# 3. run cnvnator and use copy number to annotate interspersed duplication candidates 
# (gender file is required)
# ----
# Extract copy numbers from CNVnator
mkdir -p cnvnator/log

cat gtex.lumpy.gs.svscore.low_conf.txt | cut -f 1 > dup_region.txt
echo "exit" >> dup_region.txt

source /gsc/pkg/root/root/bin/thisroot.sh
#GTEX root file for cnvnator
while read SAMPLE BAM 
do
bomb -m 10 -J $SAMPLE.cnvnator -o cnvnator/log/$SAMPLE.cnvnator.log -e cnvnator/log/$SAMPLE.cnvnator.log \
"cat dup_region.txt | \
/gscmnt/gc2719/halllab/src/speedseq/bin/cnvnator-multi \
-root $BAM \
-genotype 100 > cnvnator/$SAMPLE.txt"
done < bam.list


# run cn_evidence.py to annotate candidate callsets (whether the interDUP candidate has copy number evidence support or not)
# Test for copy number evidence (svclassfier) 
python ../scripts/cn_evidence.py \
-i gtex.lumpy.gs.svscore.low_conf.vcf \
-v gtex.lumpy.gs.svscore.low_conf.txt \
-c cnvnator \
-o gtex.lumpy.gs.svscore.low_conf.interdup.cn.vcf \
-g gender.txt \
-m 'large_sample'

#-i gtex.lumpy.gs.svscore.low_conf.vcf \  Input vcf
#-v gtex.lumpy.interdup.50.txt \	Interspersed duplication candidates file (output from inter_dup.py)
#-c cnvnator \	Path to the directory of cnvnator outputs
#-o gtex.lumpy.gs.svscore.low_conf.interdup.cn.vcf \	Output vcf
#-g gender.txt \	Gender file for samples
#-m 'large_sample'	Reclassification method, options similar to sv_classifier.py (large_sample=regression based classifier)
# MORE OPTIONS CAN BE SEEN BY `python cn_evidence.py -h`

# Interspersed duplications verified by cnvnator evidence are annotated in the output vcf file
# Following information is added to the original BND pair that re-classified as interDUP:
# [INFO]BNDTYPE (iDUP_$TYPE, $TYPE=DIRECT/INVERTED)
# [INFO]DUPPOS (coordinates of duplication region)
# [INFO]MATEID (the other BND in the same pair was added)
# [FORMAT$GT] Updated as the GT for interspersed duplication
# [FORMAT$AB] Updated as the AB for interspersed duplication

# example: 
#1	3274113	LUMPY_BND_93596_1	N	N[1:3276458[	4404.90	PASS	
# SVTYPE=BND;STRANDS=+-:271;IMPRECISE;CIPOS=-2,1;CIEND=-1,1;CIPOS95=2,2;CIEND95=1,1;
# MATEID=LUMPY_BND_93596_2,LUMPY_BND_176028_1,LUMPY_BND_176028_2;   ## LUMPY_BND_176028 is the other BND in the pair
# EVENT=LUMPY_BND_93596;SU=506;PE=224;SR=282;ALG=PROD;AF=0.04082;EXVAR=0.0162;NSAMP=12;MSQ=367.08;
# BNDTYPE=iDUP_DIRECT; ## This BND is reclassified as iDUP (interspersed duplication), and the duplicated region was directed inserted (rather than being inverted and then inserted)
# DUPPOS=1:3276456-3276585;  ## The coordinates of duplication region on the reference genome
# INSPOS=1:3274121 ## The coordinates of insertion on the reference genome

