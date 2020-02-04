#!/bin/sh
VAR_ID=$1
TRAIT=$2
echo $VAR_ID
echo $TRAIT

ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/

INVCF=${ROOT}/2-Association_test/data/
PHENO=${INVCF}/general/qt_finnseq_20190620.renormalized.wgs.5k_ids.power80.ped

TRAIT_LIST=${ROOT}/inputs/pheno/large_ss_combined_pheno_201610.list


# 2. Make plots

## boxplot & scatter plot of geno vs pheno
G_vs_P_script=${ROOT}/3-Candidate_analysis/scripts/2.1-G_vs_P.per_var.R

mkdir -p ${VAR_ID}/plots/${TRAIT}
Rscript ${G_vs_P_script} --geno=${VAR_ID}/data/genotype.t.list \
--pheno=${PHENO} --trait=${TRAIT} --out=${VAR_ID}/plots/${TRAIT} \
--assoc=${VAR_ID}/data/assoc.results.all_traits.sorted.txt 

## Phenomewide Manhattan plot
Trait_group=${ROOT}/general_info/large_ss_combined_pheno_201610.csv
Phewide_Mplot_script=${ROOT}/3-Candidate_analysis/scripts/2.2-Phewide_Mplot.per_ver.R

Rscript ${Phewide_Mplot_script} --group=${Trait_group} \
--assoc=${VAR_ID}/data/assoc.results.all_traits.sorted.txt \
--out=${VAR_ID}/plots/

