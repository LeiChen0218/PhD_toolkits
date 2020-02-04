#!/usr/bin/env python

import argparse, sys
#import pandas as pd 
# import math, time, re
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
    parser.add_argument('-a', '--alignment', required=True, type=argparse.FileType('r'), default=None, help=' alignment output (pair format)')
    # parse the arguments
    args = parser.parse_args()

    return args

def convert_align_format(align):
    "Covert alignment to table format"
    id1=None
    id2=None
    #df = pd.DataFrame(columns= ['alignment_pos','str1','str2','matched','pos1','pos2'])

    i = 0  # trace the trimer
    alignment_pos=1
    for line in align:
        if(line.startswith("#") or len(line.rstrip()) == 0):
            continue
        else:
            if i == 1: # skip the middle line
                i += 1
            else: # not middle line
                l = line.rstrip().split()
                if i == 0: # first sequence
                    if id1 is None:
                        id1 = l[0]
                    start_pos1= int(l[1])
                    seq1=l[2]
                    i += 1

                if i == 2: # second sequence
                    if id2 is None:
                        id2 = l[0]
                        header=['alignment_pos',id1,id2,'matched','pos1','pos2']
                        sys.stdout.write('\t'.join(header)+'\n')
                    start_pos2= int(l[1])
                    seq2=l[2]
                    i = 0 # start point for next group

                    # record the alignmet results
                    if len(seq1) != len(seq2):
                        sys.stderr.write('Aligned sequneces are not equal in length:')
                        sys.stderr.write(start_pos1+'/'+start_pos2+'\n')
                        continue
                    else:
                        for j in range(0,len(seq1)):
                            match = (seq1[j] == seq2[j])
                            # df.append({'alignment_pos': alignment_pos,
                            #     'str1': seq1[j],
                            #     'str2': seq2[j],
                            #     'matched': match,
                            #     'pos1': start_pos1,
                            #     'pos2': start_pos2})
                            row_list = [alignment_pos, seq1[j], seq2[j], match, start_pos1,start_pos2]
                            row_str = '\t'.join([str(x) for x in row_list])
                            sys.stdout.write(row_str+'\n')

                            if seq1[j] != '-':
                                start_pos1 += 1
                            if seq2[j] != '-':
                                start_pos2 += 1
                            alignment_pos += 1
    sys.stdout.close()
    sys.stderr.close()


# --------------------------------------
# main function

def main():
    # parse the command line args
    # parse the command line args
    args = get_args()
    convert_align_format(args.alignment)
    

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 