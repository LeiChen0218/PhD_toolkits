#!/usr/bin/env python

import argparse, sys
from scipy import stats
import gzip
from argparse import RawTextHelpFormatter

__author__ = "Lei Chen (leichen@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2019-06-13 15:50 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcf_allele_freq.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Add allele frequency information to a Genome STRiP VCF file\n\
    Allele freq. is defined as the fraction of samples that deviate from the mode copy number")
    parser.add_argument('-v', '--vcf', metavar='VCF', required=True, type=str, default=None, help='vcf with the second GT Field as CN (.gz)')
    # parse the arguments
    args = parser.parse_args()

    return args

def get_svtype(info):
    info_list=info.split(';')
    cnvtype="NA"
    for i_str in info_list:
        t = i_str.split('=')
        if t[0] == "GSCNCATEGORY":
            cnvtype = t[1]
    return cnvtype

class Variant(object):
    def __init__(self, var_list):
        self.var_id = var_list[2]
        self.svtype = get_svtype(var_list[7])
        self.cns = [ int(gt_str.split(':')[1]) for gt_str in var_list[9:len(var_list)] ]
        self.mode_cn = self.get_mode_cn()

    def get_mode_cn(self):
        return max(stats.mode(self.cns)[0])

    def carrier_list(self, ref_cn):
        return [ cn != ref_cn for cn in self.cns]

class Count(object):
    def __init__(self,sample):
        self.s=sample
        self.c_dup=0
        self.c_del=0
        self.c_mixed=0

    def add_count(self,svtype):
        if svtype == "DUP":
            self.c_dup += 1
        if svtype == "DEL":
            self.c_del += 1
        if svtype == "MIXED":
            self.c_mixed += 1

# primary function
def cal_freq_and_counts(vcf_file):

    vcount_by_sample_r2 = []
    vcount_by_sample_rm = []
    sample_list = []

    temp = 1
    # read input VCF
    for line in vcf_file:
        if line.startswith('##'):
            continue

        v = line.rstrip().split('\t')
        if line.startswith('#CHROM'):
            sample_list = v[9:len(v)]
            for s in sample_list:
                c1 = Count(s)
                vcount_by_sample_rm.append(c1)
                c2 = Count(s)
                vcount_by_sample_r2.append(c2)
            continue

        var = Variant(v)
        if 'X' in v[0] or 'Y' in v[0]: # for QC purpose, only count autosome SVs
            continue
        if var.svtype == "NA":
            continue

        ref_carrier = var.carrier_list(2)
        var_carrier = var.carrier_list(var.mode_cn)

        for i in range(len(sample_list)): 
            if ref_carrier[i]: 
                vcount_by_sample_r2[i].add_count(var.svtype)
            if var_carrier[i]:
                vcount_by_sample_rm[i].add_count(var.svtype)

    #output variant count by samples to stdout
    sys.stdout.write('sample\tDUP_r2\tDEL_r2\tMIXED_r2\tDUP_rm\tDEL_rm\tMIXED_rm\n')
    for s in range(len(sample_list)):
        #print vcount_by_sample_r2[s]
        #continue
        sys.stdout.write('\t'.join([sample_list[s],
                                str(vcount_by_sample_r2[s].c_dup),
                                str(vcount_by_sample_r2[s].c_del),
                                str(vcount_by_sample_r2[s].c_mixed),
                                str(vcount_by_sample_rm[s].c_dup),
                                str(vcount_by_sample_rm[s].c_del),
                                str(vcount_by_sample_rm[s].c_mixed)]))
        sys.stdout.write('\n')

    sys.stdout.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.vcf.endswith('.gz'):
        vcf = gzip.open(args.vcf, 'rb')
    else:
        vcf = open(args.vcf, 'r')

    # call primary function
    cal_freq_and_counts(vcf)

    # close the files
    vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 