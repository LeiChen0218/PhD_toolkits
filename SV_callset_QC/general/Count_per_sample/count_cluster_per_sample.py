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

class Sample(object):
    
    def __init__(self, sid):
        self.id = sid
        self.count_nonref = 0
        self.count_nonvar = 0
        self.cluster_nonref = []
        self.cluster_nonvar = []


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
var_gt_corr.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: correlate variants and genotypes")
    parser.add_argument('-v', '--vcf', metavar='VCF', required=True, type=str, default=None, help='VCF')
    parser.add_argument('-c', '--cluster', dest='cluster', metavar='FILE', required=True, type=argparse.FileType('r'), help='cluster bed file with the last column recording ')
    parser.add_argument('-g', '--gt', dest='gt', type=str, default='GT' , help='GT field to use ')
    
    #parse the arguments
    args = parser.parse_args()
    return args

def get_cn(gt_list):
    "extract CN from the VCF line"
    return [int(x.split(':')[1]) for x in gt_list[9:len(gt_list)] ]

def gt_to_cn(gt_str):
    "return GT to CN"
    if '.' in gt_str or gt_str == "0/0":
        return 2
    elif gt_str == "1/1":
        return 0
    else:
        return 1

def get_cn_from_gt(gt_list):
    "extract CN from GT"
    return [gt_to_cn(x.split(':')[0]) for x in gt_list[9:len(gt_list)] ]


def  count_cluster_per_sample(vcf, cluster_dic,gt_field):
    sys.stdout.write('SAMPLE\tCOUNT_NONREF\tCOUNT_NONVAR\n')

    sample_count = []

    for line in vcf:
        if line.startswith('##'):
            continue
        l = line.rstrip().split('\t')
        if l[0].startswith('#'):
            sample_count = [ Sample(s) for s in l[9:len(l)]]

        if l[2] in cluster_dic.keys():
            if gt_field == "CN":
                cn_list = get_cn(l)
            else:
                cn_list = get_cn_from_gt(l)
            mode_cn = stats.mode(cn_list).mode[0]
            cluster = cluster_dic[l[2]]
            for i in range(len(cn_list)):
                if cn_list[i] != 2 and cluster not in sample_count[i].cluster_nonref :
                    sample_count[i].cluster_nonref.append(cluster)
                    sample_count[i].count_nonref += 1
                if cn_list[i] != mode_cn and cluster not in sample_count[i].cluster_nonvar :
                    sample_count[i].cluster_nonvar.append(cluster)
                    sample_count[i].count_nonvar += 1

    for samp in sample_count:
        sys.stdout.write('\t'.join( [samp.id, str(samp.count_nonref), str(samp.count_nonvar)] ))
        sys.stdout.write('\n')

    sys.stdout.close()


# --------------------------------------
# main function

def main():
    # parse the command line args
    # parse the command line args
    args = get_args()

    # get list of variants to examine
    cluster_dic={}
    if args.cluster:
        for line in args.cluster:
            v = line.rstrip().split('\t')
            cluster_dic[v[3]] = v[-1]
        args.cluster.close()
        sys.stderr.write('finished loading the cluster list\n')

    # open the VCF files
    if args.vcf.endswith('.gz'):
        vcf = gzip.open(args.vcf, 'rb')
    else:
        vcf = open(args.vcf, 'r')

    count_cluster_per_sample(vcf, cluster_dic, args.gt)

    sys.stderr.close()
    

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 