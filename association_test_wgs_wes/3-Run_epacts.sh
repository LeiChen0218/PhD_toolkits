# 3-Run_epacts.sh
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
DIR=$ROOT/2-Association_test

cd $DIR/data

LUMPY=$DIR/data/lumpy/lumpy.highQ.mac10.autosome.dechr.vcf.gz
GS=$DIR/data/gs/gs.highQ.mac10.autosome.dechr.vcf.gz
CNVNATOR=$DIR/data/cnvnator/cnvnator.highQ.mac10.autosome.dechr.ref.vcf.gz

PHENO=$DIR/data/general/qt_finnseq_20190620.renormalized.wgs.5k_ids.power80.ped

head -n 1 general/qt_finnseq_20190620.renormalized.wgs.5k_ids.power80.ped \
| sed 's/\t/\n/g' | grep "_rn"  --color=never >  general/renormalized.trait.list

trait_list=$DIR/data/general/renormalized.trait.list
KIN_MAT=/gscmnt/gc2802/halllab/lganel/FinMetSeq5k/kinship.kinf
EPACTS=/gscmnt/gc2719/halllab/users/lganel/src/EPACTS-3.2.6/scripts/epacts
export LD_LIBRARY_PATH=/gscmnt/gc2719/halllab/src/gcc-4.9.2/lib64:$LD_LIBRARY_PATH

while read trait
do
mkdir -p gs/epacts_p80qt/${trait}/logs
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo gs/epacts_p80qt/${trait}/logs/finmetseq.assoc.${trait}.log \
"${EPACTS} single \
--vcf ${GS} \
--ped ${PHENO} \
--field CN \
--min-maf -1000000 \
--min-mac -1000000 \
--kin ${KIN_MAT} \
--pheno ${trait} \
--test q.emmax \
--out gs/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc \
--run 5"
done < $trait_list

while read trait
do
	ll  gs/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc.epacts.gz
done < $trait_list 

while read trait
do
mkdir -p cnvnator/epacts_p80qt/${trait}/logs
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo cnvnator/epacts_p80qt/${trait}/logs/finmetseq.assoc.${trait}.log \
"${EPACTS} single \
--vcf ${CNVNATOR} \
--ped ${PHENO} \
--field CNF \
--min-maf -1000000 \
--min-mac -1000000 \
--kin ${KIN_MAT} \
--pheno ${trait} \
--test q.emmax \
--out cnvnator/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc \
--run 5"
done < $trait_list

while read trait
do
	ll  cnvnator/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc.epacts.gz
done < $trait_list 


while read trait
do
mkdir -p lumpy/epacts_p80qt/${trait}/logs
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo lumpy/epacts_p80qt/${trait}/logs/finmetseq.assoc.${trait}.log \
"${EPACTS} single \
--vcf ${LUMPY} \
--ped ${PHENO} \
--field AB \
--min-maf -1000000 \
--min-mac -1000000 \
--kin ${KIN_MAT} \
--pheno ${trait} \
--test q.emmax \
--out lumpy/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc \
--run 5"
done < $trait_list

while read trait
do
	ll  lumpy/epacts_p80qt/${trait}/finmetseq.5k.$trait.emmax.kin.lm.assoc.epacts.gz
done < $trait_list 

#####
# count # carriers in each phenotype 
count_scr=$DIR/scripts/general/3a-count_af_pheno.py
cd general
#trait="Alb_combined_rn"
while read trait
do
mkdir -p ac_by_pheno/${trait}/logs
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo ac_by_pheno/${trait}/logs/gs.${trait}.log \
bash -c \
"python $count_scr -v ${GS} -p ${PHENO} -g CN -t $trait > ac_by_pheno/${trait}/gs.mac.${trait}.txt"
done < $trait_list

while read trait
do
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo ac_by_pheno/${trait}/logs/lumpy.${trait}.log \
bash -c \
"python $count_scr -v ${LUMPY} -p ${PHENO} -g GT -t $trait > ac_by_pheno/${trait}/lumpy.mac.${trait}.txt"
done < $trait_list

while read trait
do
bsub -n 4 -M 10000000 \
-g /lchen/SVassoc \
-a 'docker{registry.gsc.wustl.edu/genome/genome_perl_environment}' \
-q ccdg \
-oo ac_by_pheno/${trait}/logs/cnvnator.${trait}.log \
bash -c \
"python $count_scr -v ${CNVNATOR} -p ${PHENO} -g CN -t $trait > ac_by_pheno/${trait}/cnvnator.mac.${trait}.txt"
done < $trait_list


#while read trait ; do
#	wc -l ac_by_pheno/${trait}/*txt
#done < $trait_list | awk '{print $1}' | sort | uniq -c

while read trait ; do
	cat ac_by_pheno/${trait}/*.mac.*.txt | grep -v "ID" | awk '$2>9' > ac_by_pheno/${trait}/var.mac10.${trait}.txt
done < $trait_list 

while read trait ; do
zcat lumpy/epacts_p80qt/${trait}/finmetseq.*.assoc.epacts.gz | cut -f 4 | sed 's/_/\t/2' \
| cut -f 2 | sed 's/MARKER_ID/ID/' | paste - <(zcat lumpy/epacts_p80qt/$trait/finmetseq.*epacts.gz) \
| sed 's/ /\t/g' | zjoin -a stdin -b general/ac_by_pheno/${trait}/var.mac10.${trait}.txt -wa -p "ID" \
> lumpy/epacts_p80qt/${trait}/finmetseq.${trait}.mac10.assoc.epacts
done < $trait_list 

while read trait ; do
zcat gs/epacts_p80qt/${trait}/finmetseq.*.assoc.epacts.gz | cut -f 4 | sed 's/_/\t/2' \
| cut -f 2 | sed 's/MARKER_ID/ID/' | paste - <(zcat gs/epacts_p80qt/$trait/finmetseq.*epacts.gz) \
| sed 's/ /\t/g' | zjoin -a stdin -b general/ac_by_pheno/${trait}/var.mac10.${trait}.txt -wa -p "ID" \
> gs/epacts_p80qt/${trait}/finmetseq.${trait}.mac10.assoc.epacts
done < $trait_list

while read trait ; do
zcat cnvnator/epacts_p80qt/${trait}/finmetseq.*.assoc.epacts.gz | cut -f 4 | sed 's/_/\t/2' \
| cut -f 2 | sed 's/MARKER_ID/ID/' | paste - <(zcat cnvnator/epacts_p80qt/$trait/finmetseq.*epacts.gz) \
| sed 's/ /\t/g' | zjoin -a stdin -b general/ac_by_pheno/${trait}/var.mac10.${trait}.txt -wa -p "ID" \
> cnvnator/epacts_p80qt/${trait}/finmetseq.${trait}.mac10.assoc.epacts
done < $trait_list 






