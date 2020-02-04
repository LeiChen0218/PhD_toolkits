#!/bin/bash

set -ueo pipefail

VCF=$1
VAWK=/gscmnt/gc2719/halllab/bin/vawk

# NOTE - we're removing secondary variants from this file so we don't double count them in plink
zcat $VCF | \
	$VAWK --header '{ \
	split(I$STRANDS,x,","); \
	split(x[1],y,":"); \
	split(x[2],z,":"); \
    if (I$SVTYPE=="DEL" || I$SVTYPE=="DUP" || I$SVTYPE=="MEI"){ \
	$7="PASS"; print $0; \
	}  else if ( I$SVTYPE=="INV" && I$MSQ>=100 && (I$SR/I$SU)>=0.1 && (I$PE/I$SU)>=0.1 && (y[2]/I$SU)>0.1 && (z[2]/I$SU)>0.1){ \
	$7="PASS"; print $0; \
    } else if ( I$SVTYPE=="BND" && I$MSQ>=250){ \
	$7="PASS"; print $0; \
	} else { \
	$7="LOW"; print $0; \
	} \
}' 
