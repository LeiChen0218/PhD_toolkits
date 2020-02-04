#!/bin/bash
x=$1
sdir=$2
region_file=$3
out_dir=$4
CNVNATOR_MULTI=/opt/ccdg/cnvnator-0.3.3/bin/cnvnator

echo "SNP,Chromosome,PhysicalPosition,$x" | sed 's/,/\t/g' > ${out_dir}/cn/$x.cn
cat $region_file |\
${CNVNATOR_MULTI} -root $sdir/*cnvnator.out/*.cram.hist.root \
-genotype 100 |\
awk '{OFS="\t"; if($1 !~ /male/){print $2,$4} }' |\
sed s/[:-]/"\t"/g |\
awk '{OFS="\t"; print $1":"$2"-"$3,$1,$2,$4}' >> ${out_dir}/cn/$x.cn
