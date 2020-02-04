#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from argparse import RawTextHelpFormatter

__author__ = "Lei Chen (leichen@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2018-05-14e 14:53 $"

# --------------------------------------
# define functions


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
var_gt_corr.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: correlate variants and genotypes")
    parser.add_argument('-v', '--vcf', metavar='VCF', required=True, type=str, default=None, help='VCF')
    parser.add_argument('-p', '--ped', dest='ped', metavar='FILE', required=True, type=argparse.FileType('r'), help='Phenotype file in .ped format')
    parser.add_argument('-g', '--gt', dest='gt', type=str, default='GT' , help='GT field to use ')
    parser.add_argument('-t', '--trait', dest='trait', type=str, required=True, help='trait name')
    #parse the arguments
    args = parser.parse_args()
    return args

def get_cn(gt_list, gt_index):
    "extract CN from the VCF line"
    l = []
    #print gt_index
    for i in gt_index:
        gt = gt_list[i].split(':')
        gt_cn = int(gt[1])
        l.append(gt_cn)
    return l

def gt_to_cn(gt_str):
    "return GT to CN"
    if '.' in gt_str:
        return np.nan
    if gt_str == "0/0":
        return 2
    elif gt_str == "1/1":
        return 0
    else:
        return 1

def get_cn_from_gt(gt_list, gt_index):
    "extract CN from GT"
    l = []
    for i in gt_index:
        gt = gt_list[i].split(':')
        gt_cn = gt_to_cn(gt[0])
        l.append(gt_cn)
    return l

def process_ped(ped, trait):
    sample_list = []
    #trait_index = 0
    for line in ped:
        l = line.rstrip().split('\t')
        if l[0] == "#FAM_ID":
            trait_index = l.index(trait)
            #print trait_index
        else:
            if trait_index == 0:
                sys.stderr.write("Can not find query trait in input ped file")
                return 
            if l[trait_index] != "NA":
                sample_list.append(l[0])
    #print len(sample_list)
    return sample_list


def count_ac_per_pheno(vcf, sample_list, gt_field):
    sys.stdout.write('ID\tAC\tAF\tN\n')
    sample_index = []

    for line in vcf:
        if line.startswith('##'):
            continue

        l = line.rstrip().split('\t')
        if l[0].startswith('#'):
            for i in range(len(l)):
                if l[i] in sample_list:
                    sample_index.append(i)
            continue

        if gt_field == "CN":
            cn_list = get_cn(l, sample_index)
            #print cn_list
        else:
            cn_list = get_cn_from_gt(l, sample_index)

        mode_cn = stats.mode(cn_list).mode[0]
        non_missing = sum([x is not np.nan for x in cn_list])
        mode_count =  sum([x == mode_cn for x in cn_list])
        mac = non_missing - mode_count
        maf = float(mac)/non_missing

        sys.stdout.write('\t'.join([ l[2],str(mac),str(maf), str(non_missing) ]))
        sys.stdout.write('\n')

    sys.stdout.close()


# --------------------------------------
# main function

def main():
    # parse the command line args
    # parse the command line args
    args = get_args()

    sample_list = []
    if not (args.ped and args.trait):
        sys.exit('Missing phenotype file or trait')
    else:
        sample_list = process_ped(args.ped,args.trait)

    # open the VCF files
    if args.vcf.endswith('.gz'):
        vcf = gzip.open(args.vcf, 'rb')
    else:
        vcf = open(args.vcf, 'r')

    count_ac_per_pheno(vcf, sample_list, args.gt)

    sys.stderr.close()
    

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 