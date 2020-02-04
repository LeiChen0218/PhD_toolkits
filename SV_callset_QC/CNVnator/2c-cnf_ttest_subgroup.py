#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import numpy as np
from scipy import stats
from decimal import Decimal

__author__ = "Lei Chen (leichen@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2019-06-22 09:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcf_allele_freq.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: calculate allele frequency information in two subset of samples from a Genome STRiP VCF file\n\
    Allele freq. is defined as the fraction of samples that deviate from the mode copy number")
    parser.add_argument('-v',metavar='vcf', dest='input_vcf', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')
    parser.add_argument('-s','--subset',metavar='<FILE>',dest='subset',type=argparse.FileType('r'),default=None, help='subset of samples')
    parser.add_argument('-g','--gt_field',metavar='<FILE>',dest='gt_field',type=str,default="CN", help='CN field used for t-test (default: CN)')
        # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = sys.stdin
    # send back the user input
    return args


def get_gt(format_str, cn_index):
    l = format_str.split(':')
    gt = l[cn_index]
    if gt == '.':
        return np.nan
    return float(gt)


# primary function
def cal_freq_and_counts(vcf_file,subset,gt_field):
    
    output = sys.stdout
    output.write("SV_id\tt\tp-value\n")

    sample_list = []
    # read input VCF
    for line in vcf_file:
        if line.startswith('##'):
             continue
        v = line.rstrip().split('\t')
        if line.startswith('#CHROM'):
            sample_list = v[9:len(v)]    
            continue

        format=v[8].split(':')
        if gt_field not in format:
            sys.stderr.write("Warning: input gt field does not exist, use default field CN")
            gt_field = "CN"
        gt_index = format.index(gt_field)
        #print gt_index

        # 1- carrier, 0-non carrier
        gt_list = []
        subset_gt_list = []
        nonsubset_gt_list = []

        gt_cols = v[9:len(v)]
        for i in range(len(sample_list)):
            cn = get_gt(gt_cols[i], gt_index)
            gt_list.append(cn)
            s = sample_list[i]
            if s in subset:
                subset_gt_list.append(cn)
            else:
                nonsubset_gt_list.append(cn)
            
        t, prob = stats.ttest_ind(subset_gt_list, nonsubset_gt_list, equal_var=False )
        
        output.write( '\t'.join([v[2],
                                  '%0.4f' % t,
                                  '%0.4E' % Decimal(prob)]))
        output.write('\n')

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    subset=[]
    for line in args.subset:
        subset.append(line.rstrip())

    # call primary function
    cal_freq_and_counts(args.input_vcf,subset,args.gt_field)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 