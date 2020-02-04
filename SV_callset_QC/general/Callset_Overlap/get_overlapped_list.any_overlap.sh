#!/bin/sh

bed_a=$1 #output list, 4th column as the var_id 
bed_b=$2 #query bed
output=$3 #output list

bedtools=/gscmnt/gc2719/halllab/bin/bedtools

# reciprocal overlap > 0.5
${bedtools} intersect -a $bed_a -b $bed_b | cut -f 4 | sort | uniq > tmp.list

# combine 
zjoin -a $bed_a -b tmp.list -14 -r -e "No" | awk '{OFS="\t";print $4,$NF}' \
| awk '{OFS="\t";if($2=="No")print; else print $1,"Yes"}' \
> $output

rm -f tmp.list