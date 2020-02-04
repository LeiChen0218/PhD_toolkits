#!/bin/sh

bed_a=$1 #output list, 4th column as the var_id 
bed_b=$2 #query bed
output=$3 #output list

bedtools=/gscmnt/gc2719/halllab/bin/bedtools

# reciprocal overlap > 0.5
${bedtools} intersect -a $bed_a -b $bed_b -r -f 0.5 | cut -f 4 | sort | uniq > tmp.r05.list

# small variant inside big ones in the query
${bedtools} intersect -a $bed_a -b $bed_b -f 1 | cut -f 4 | sort | uniq > tmp.f08.list

# combine 
cat tmp.r05.list tmp.f08.list | sort | uniq \
| zjoin -a $bed_a -b stdin -14 -r -e "No" | awk '{print $4,$NF}' \
| awk '{OFS="\t";if($2=="No")print; else print $1,"Yes"}' \
> $output

rm -f tmp.r05.list
rm -f tmp.f08.list