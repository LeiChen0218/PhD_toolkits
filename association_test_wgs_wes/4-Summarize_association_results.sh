# 4-Summarize_association_results.sh (re-ran for MAC10)
ROOT=/gscmnt/gc2802/halllab/lchen/finmetseq_final_20190604/
DIR=$ROOT/2-Association_test

##
# get the epacts results from cluster and make plots on local
#

# get candidates

# wes replicate : 0.001
# genome-wide : 1.58E-06

cd $DIR/data

trait_list=$DIR/data/general/renormalized.trait.list
#trait="Alb_combined_rn"
printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\n" > lumpy/lumpy.candidate.p3.txt
while read trait ; do
awk '$12<0.001' lumpy/epacts_p80qt/${trait}/finmetseq.*.mac10.*.epacts | cut -f 1,2,3,12,13 \
| zjoin -a stdin -b general/ac_by_pheno/${trait}/lumpy.*.txt | cut -f 6 --complement \
| awk -v t=$trait '{OFS="\t";print t,$0}'
done < $trait_list >> lumpy/lumpy.candidate.p3.txt

printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\n" > cnvnator/cnvnator.candidate.p3.txt
while read trait ; do
awk '$12<0.001' cnvnator/epacts_p80qt/${trait}/finmetseq.*.mac10.*.epacts | cut -f 1,2,3,12,13 \
| zjoin -a stdin -b general/ac_by_pheno/${trait}/cnvnator.*.txt | cut -f 6 --complement \
| awk -v t=$trait '{OFS="\t";print t,$0}'
done < $trait_list >> cnvnator/cnvnator.candidate.p3.txt

printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\n" > gs/gs.candidate.p3.txt
while read trait ; do
awk '$12<0.001' gs/epacts_p80qt/${trait}/finmetseq.*.mac10.*.epacts | cut -f 1,2,3,12,13 \
| zjoin -a stdin -b general/ac_by_pheno/${trait}/gs.*.txt | cut -f 6 --complement \
| awk -v t=$trait '{OFS="\t";print t,$0}'
done < $trait_list >> gs/gs.candidate.p3.txt

##
#printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\tMETHOD\n" > general/candidate.1E6.txt
#awk '{OFS="\t";if($5<0.00000158) print $0,"LUMPY"}' lumpy/lumpy.candidate.p3.txt >> general/candidate.1E6.txt
#awk '{OFS="\t";if($5<0.00000158) print $0,"GS"}' gs/gs.candidate.p3.txt >> general/candidate.1E6.txt
#awk '{OFS="\t";if($5<0.00000158) print $0,"CNVNATOR"}' cnvnator/cnvnator.candidate.p3.txt >> general/candidate.1E6.txt

# MAC10: 26,495.3 independent variables
# new genome-wide-sig thres: 
# 
printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\tMETHOD\n" > general/candidate.2E6.txt
awk '{OFS="\t";if($5<0.00000189) print $0,"LUMPY"}' lumpy/lumpy.candidate.p3.txt >> general/candidate.2E6.txt
awk '{OFS="\t";if($5<0.00000189) print $0,"GS"}' gs/gs.candidate.p3.txt >> general/candidate.2E6.txt
awk '{OFS="\t";if($5<0.00000189) print $0,"CNVNATOR"}' cnvnator/cnvnator.candidate.p3.txt >> general/candidate.2E6.txt


## get sub-threshold variants 
printf "TRAIT\tID\tCHR\tPOS\tP\tBETA\tAC\tAF\tN\tMETHOD\n" > general/candidate.1E5.txt
awk '{OFS="\t";if($5<0.00001) print $0,"LUMPY"}' lumpy/lumpy.candidate.p3.txt >> general/candidate.1E5.txt
awk '{OFS="\t";if($5<0.00001) print $0,"GS"}' gs/gs.candidate.p3.txt >> general/candidate.1E5.txt
awk '{OFS="\t";if($5<0.00001) print $0,"CNVNATOR"}' cnvnator/cnvnator.candidate.p3.txt >> general/candidate.1E5.txt


