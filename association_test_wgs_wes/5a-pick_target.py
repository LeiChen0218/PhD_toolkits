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
__date__ = "$Date: 2019-07-30e 14:53 $"

class exonicCNV(object):

    def __init__(self, cnv, exons):
        self.cnv = cnv
        self.exons = exons.split(',')
        self.cnv_dosage = None
        #self.exon_r2 = {}
        #self.exon_r2 = self.init_exon_r2(self.exons)

    def init_cnv_dosage(self, dosage, sample, samp_dic):
        sample_oid = []
        for s in sample:
            sample_oid.append(samp_dic[s])
        self.cnv_dosage = pd.Series(dosage, index=sample_oid, name="cnv")
        return

    def cal_r2(self, exon_dosage):

        if self.cnv_dosage is None:
            sys.stderr.write('Warning: dosage is missing for CNV: '+self.cnv)
            sys.stderr.write('\n')
            return -1
        df = pd.concat([self.cnv_dosage,exon_dosage],axis=1,join='inner').dropna()
        #print len(df)
        if len(set(df["cnv"])) == 1:
            return -1 
        slope, intercept, r_value, p_value, std_err = stats.linregress(df["cnv"],df["exon"])

        return r_value**2

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
var_gt_corr.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: correlate variants and genotypes")
    parser.add_argument('-l', '--list', metavar='FILE', required=True, type=argparse.FileType('r'), default=None, help='list of CNV and its exon candidates separated by comma')
    parser.add_argument('-v', '--sv', metavar='FILE', required=True, type=argparse.FileType('r'), default=None, help='SV dosage file')
    parser.add_argument('-e', '--exon', dest='exon', metavar='FILE', required=True, type=argparse.FileType('r'), help='Exon dosage file')
    parser.add_argument('-s', '--samp_dic', metavar='FILE', dest='samp_dic', type=argparse.FileType('r'), default=None, required=True, help='sample id dictionary')
    parser.add_argument('-o', '--out', metavar='FILE', dest='out', type=argparse.FileType('w'), default=None, required=True, help='output file')

    # parse the arguments
    args = parser.parse_args()

    return args

def read_dosage(dosage_str):
    if dosage_str == '.':
        return np.nan
    else:
        return float(dosage_str)


def calculate_r2_sv_exons(exons, SVs, out):
    in_header = True
    samp_list = []
    exon_dosage_dic = {}

    for line in exons:
        l = line.rstrip().split('\t')
        if in_header:
            samp_list = l[1:len(l)]
            in_header = False
        else:
            dosage = [read_dosage(x) for x in l[1:len(l)]]
            exon_dosage = pd.Series(dosage, index=samp_list, name="exon")
            exon_dosage_dic[l[0]] = exon_dosage
    sys.stderr.write('finished loading the exon dosage\n')

    count = 1
    cnv_count = len(SVs.keys())
    for sv in SVs.values():
        sys.stderr.write('Progress: '+str(count)+'/'+str(cnv_count))
        sys.stderr.write('\n')
        count += 1
        for exon in sv.exons:
            if exon not in exon_dosage_dic.keys():
                continue 
            exon_dosage = exon_dosage_dic[exon]
            r2 = sv.cal_r2(exon_dosage)
            #print '\t'.join([sv.cnv, exon, "%.3f" % r2])
            out.write( '\t'.join([sv.cnv, exon, "%.3f" % r2]))
            out.write('\n')

    out.close()


# --------------------------------------
# main function

def main():
    # parse the command line args
    # parse the command line args
    args = get_args()

    # SV -> exons
    SVs = {}
    if args.list is not None:
        for line in args.list:
            l = line.rstrip().split('\t')
            sv = exonicCNV(l[0],l[1])
            SVs[l[0]] = sv
            #sys.stderr.write(sv.cnv)
            #sys.stderr.write('\n')
        args.list.close()
        sys.stderr.write('finished loading the variants list\n')
    else:
        parser.print_help()
        exit(1)

    # get sample dictionary
    samp_dic = {}
    if args.samp_dic is not None:
        for line in args.samp_dic:
            v = line.rstrip().split('\t')
            samp_dic[v[0]] = v[1] # nid -> oid
        args.samp_dic.close()
        sys.stderr.write('finished loading the sample dics\n')
    else:
        parser.print_help()
        exit(1)

    # initiate cnv dosage info
    samp_list = []
    if args.sv is not None :
        in_header = True
        for line in args.sv:
            l = line.rstrip().split('\t')
            if in_header:
                samp_list = l[1:len(l)]
                in_header = False
            else:
                if l[0] in SVs.keys():
                    dosage = [read_dosage(x) for x in l[1:len(l)]]
                    SVs[l[0]].init_cnv_dosage(dosage, samp_list, samp_dic)
        sys.stderr.write('finished loading the cnv dosage\n')
    else:
        parser.print_help()
        exit(1)

    if args.exon is not None:
        calculate_r2_sv_exons(args.exon, SVs, args.out)
    else:
        parser.print_help()
        exit(1) 
    

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 