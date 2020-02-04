#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse, sys, math
from scipy import stats
import gzip

def gt_to_dosage(GT):
        "Convert genotypes to numbers for calculating correlation coefficient"
        if '/' in GT:
            if GT == '0/0':
                return 0
            elif GT == '0/1':
                return 1
            elif GT == '1/1':
                return 2
            else:
                return np.nan
        else:
            return float(GT)
 
def get_gt_list(vcf_cols):

    return [ x.split(':')[0] for x in vcf_cols ]

def is_biallelic(gt_list):
    "See if the SNP is biallelic"
    alleles=set(gt_list)
    for allele in alleles:
        if allele not in ['0/0','0/1','1/1','./.']:
            return False
    return True

def calculate_pairwise_ld(GT1,GT2):

    slope, intercept, r_value, p_value, std_err = stats.linregress(GT1,GT2)
    
    return r_value**2

def calculate_ld(neighbor_g,hit_g,output):
    ld_out = open(output,'w')
    ld_out.write("snp1\tsnp2\tdprime\trsquare\n")

    target_sample_list = []
    in_header = True
    for line in hit_g:
        l = line.rstrip().split('\t')
        if in_header:
            target_sample_list = l[1:]
            in_header = False
        else:
            target = l[0]
            dosage_list = [gt_to_dosage(x) for x in l[1:]]
            target_geno = pd.Series(dosage_list, index=target_sample_list,name=target)
    
    neighbor_sample_list = []
    for line in neighbor_g:
        if line.startswith('##'):
            continue
        else:
            l = line.rstrip().split('\t')
            if l[0].startswith("#"):
                neighbor_sample_list = l[9:]
            else:
                var_name = ':'.join(l[0:2])
                dosage_list = []
                gt_list = get_gt_list(l[9:])
                #print gt_list
                if is_biallelic(gt_list) == False:
                    continue
                dosage_list = [gt_to_dosage(x) for x in gt_list]

                neighbor_geno = pd.Series(dosage_list, index=neighbor_sample_list,name=var_name)
                
                #print target_geno           
                df = pd.concat([target_geno,neighbor_geno],axis=1,join='inner').dropna()

                r_sqr = calculate_pairwise_ld(df[target],df[var_name])
                ld_out.write(var_name+' '+target+' NA '+str(r_sqr)+'\n')

    ld_out.close()

# Arguments parsing
def description():
    return 'Association test for CNV data'

def add_argument_to_parser(parser):
    
    parser.set_defaults(entry_point=run_from_args)

def get_args():
    parser = argparse.ArgumentParser(description=description())
    parser.add_argument('-t','--target',metavar='<FILE>',dest='target',type=argparse.FileType('r'),required=True,default=None,help='Genotype file for the SV target')
    parser.add_argument('-v', '--vcf', metavar='VCF', required=True, type=str, default=None, help='flanking SNP VCF file')
    parser.add_argument('-o','--output',metavar='<FILE>',dest='output',type=str,default=None,help='output LD file')
    return parser.parse_args()
    
# ------------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()
    if args.target is None:
        parser.print_help()
        exit(1)
    # open the VCF files
    if args.vcf.endswith('.gz'):
        vcf = gzip.open(args.vcf, 'rb')
    else:
        vcf = open(args.vcf, 'r')

    # call primary function
    calculate_ld(   vcf,
                    args.target,
                    args.output,)

    args.target.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
        