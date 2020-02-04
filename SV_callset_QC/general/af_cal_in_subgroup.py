#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

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
    parser.add_argument('-g','--gt_field',metavar='<FILE>',dest='gt_field',type=str,default="GT", help='GT field used for calculating AF (default: GT)')
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

def carrier_gt(gt):
    if gt == '0/0' or '.' in gt:
        return 0
    else:
        return 1

def carrier_cn(cn):
    if cn == "2":
        return 0
    else:
        return 1

def get_gt(format_str, gt_index):
    l = format_str.split(':')
    gt = l[gt_index]
    carrier_status = 0
    if '/' in gt:
        carrier_status = carrier_gt(gt)
    else:
        carrier_status = carrier_cn(gt)
    return carrier_status


# primary function
def cal_freq_and_counts(vcf_file,subset,gt_field):
    
    output = sys.stdout
    output.write("SV_id\tcarrier_rate_subset\tcarrier_rate_nonsubset\tSV_carrier_rate\n")

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
            sys.stderr.write("Warning: input gt field does not exist, use default field GT")
            gt_field = "GT"
        gt_index = format.index(gt_field)
        #print gt_index

        # 1- carrier, 0-non carrier
        gt_list = []
        subset_gt_list = []
        nonsubset_gt_list = []

        gt_cols = v[9:len(v)]
        for i in range(len(sample_list)):
            carrier_status = get_gt(gt_cols[i], gt_index)
            gt_list.append(carrier_status)
            s = sample_list[i]
            if s in subset:
                subset_gt_list.append(carrier_status)
            else:
                nonsubset_gt_list.append(carrier_status)
            
        #print v[2],len(gt_list),len(subset_gt_list),len(nonsubset_gt_list)
        # calculate allele frequency for each CNV
        if len(gt_list) == 0:
            continue
        af = float(sum(gt_list))/len(gt_list)
        if len(subset_gt_list) == 0:
            continue
        af_sub = float(sum(subset_gt_list))/len(subset_gt_list)
        if len(nonsubset_gt_list) == 0:
            continue
        af_nonsub = float(sum(nonsubset_gt_list))/len(nonsubset_gt_list)

        output.write( '\t'.join([v[2],
                                  '%0.4f' % af_sub,
                                  '%0.4f' % af_nonsub,
                                  '%0.4f' % af]))
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