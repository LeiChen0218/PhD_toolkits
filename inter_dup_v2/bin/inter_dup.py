#! /usr/bin/env python
import sys
from optparse import OptionParser
import re
from sv_bedpe_parser import VariantList, VariantDetails
from scipy import stats

# interspersed duplication class
class interspersed_dup(object):
    
    def __init__(self,dup_type):
    	self.bnd1 = None
        self.bnd2 = None
        self.chr_dup = 0
        self.chr_ins = 0
        # estimated start and end positions of duplication region
        self.start = 0
        self.end = 0
        # estimated inserted position
        self.insert = 0
        # orientation type: normal or inverted
        if dup_type == 1:
            self.dup_type = "DIRECT"
        elif dup_type == 2:
       		self.dup_type = "INVERTED"
        else:
            self.dup_type = None
        self.stats = {"pearsonr":0,"kappa-gt":0,"kappa-binary":0}
        self.samples = []
        # genotype agreement table [[0/0&0/0,0/0&1/0],[1/0&0/0,1/1&1/1]]
        self.gt_agreement = [[0,0],[0,0]]
        # estimated genotype and allele balance for samples
        self.format = {"GT":[],"AB":[]}

    def dup_start(self,startA,startB):
        if startA < startB:
            self.start = int(startA)
        else:
            self.start = int(startB)

    def dup_end(self,endA,endB):
        if endA > endB:
            self.end = int(endA)
        else:
            self.end = int(endB)

    def dup_insert(self,startA,startB,endA,endB):
        self.insert = int((startA+startB+endA+endB)/4)

    def GT_to_number(self,GT):
        "Convert genotypes to numbers for calculating correlation coefficient"
        if GT == '0/0':
            return 0
        elif GT == '0/1':
            return 1
        elif GT == '1/1':
            return 2
        else:
            return None

    def cohen_kappa_test(self,table,test):
        "Calculate Kappa statistic for the given agreement table"
        #print table
        pe = 0
        po = 0
        total = sum(table[0])+sum(table[1])+sum(table[2])
        if test == "kappa-binary":
            pa = [0,0]
            pb = [0,0]
            pa[0] = float(sum(table[0]))/total
            pa[1] = float(sum(table[1]+table[2]))/total
            pb[0] = float(table[0][0]+table[1][0]+table[2][0])/total
            pb[1] = float(table[0][1]+table[1][1]+table[2][1]+table[0][2]+table[1][2]+table[2][2])/total
            for i in [0,1]:
                pe += pa[i]*pb[i]
                po += table[i][i]
            po += table[2][2]+table[1][2]+table[2][1]
        elif test == "kappa-gt":
            pa = [0,0,0]
            pb = [0,0,0]
            for i in [0,1,2]:
                pa[i] = float(sum(table[i]))/total
                pb[i] = float(table[0][i]+table[1][i]+table[2][i])/total
                pe += pa[i]*pb[i]
                po += table[i][i]
        po = float(po)/total
        if pe == 1:
            return "NA"
        k = (po-pe)/(1-pe)
        return k

    def parse_format(self,format1,format2):
        "compare format values from two BNDs and calculate that for interspersed duplication"
        if format1 == '.' and format2 == '.':
            return '.'
        elif format1 == '.':
            return format2
        elif format2 == '.':
            return format1
        else:
            # Genotype value
            if '/' in format1:
                GT1 = self.GT_to_number(format1)
                GT2 = self.GT_to_number(format2)
                if GT1 > GT2:
                    return format1
                else:
                    return format2
            else:
                if float(format1) > float(format2):
                    return format1
                else:
                    return format2

    def AB_recalculate(self,QA1,QA2,QR1,QR2):
        "Recalculate the allele balance from info of two BNDs"
        if QA1 and QA2 and QR1 and QR2:
            if float(QR1)+float(QR2)+float(QA1)+float(QA2) > 0:
                return (float(QA1)+float(QA2))/(float(QR1)+float(QR2)+float(QA1)+float(QA2))
        elif QA1 and QR1:
            if float(QR1)+float(QA1) > 0:
                return float(QA1)/(float(QR1)+float(QA1))
        elif QA2 and QR2:
            if float(QR2)+float(QA2) > 0: 
               return float(QA2)/(float(QR2)+float(QA2))
        else:
            return '.'
                 
    def check_GT(self,BND1,BND2,sample_list):
        "Check the genotypes agreements of two BNDs across samples (using given test)"
        self.samples = sample_list
        self.format['GT'] = range(len(sample_list))
        self.format['AB'] = range(len(sample_list))
        # kappa table :  BND2
        #               0 1 2
        #             0 A D D
        #       BND1  1 D A D
        #             2 D D A
        # Store the frequencies of agreement(A) and disagreement(D)
        kappa_table = [ [0,0,0], \
                        [0,0,0], \
                        [0,0,0]]
        i = 0
        while i < len(BND1.format['GT']):
            GT1 = self.GT_to_number(BND1.format['GT'][i])
            GT2 = self.GT_to_number(BND2.format['GT'][i])
            self.format['GT'][i] = self.parse_format(BND1.format['GT'][i],BND2.format['GT'][i])
            if GT1 != None and GT2 != None:
                kappa_table[GT1][GT2] += 1
            i += 1
        self.gt_agreement[0][0] = kappa_table[0][0]
        self.gt_agreement[0][1] = kappa_table[0][1] + kappa_table[0][2]
        self.gt_agreement[1][0] = kappa_table[1][0] + kappa_table[2][0]
        self.gt_agreement[1][1] = kappa_table[1][1] + kappa_table[2][2] + kappa_table[1][2] + kappa_table[2][1]

        #self.stats["kappa-gt"] = self.cohen_kappa_test(kappa_table,"kappa-gt")
        #self.stats["kappa-binary"] = self.cohen_kappa_test(kappa_table,"kappa-binary   ")

        # Calculate pearson R
        AB1 = []
        AB2 = []
        i = 0
        while i < len(BND1.format['AB']):
            self.format['AB'][i] = self.AB_recalculate(BND1.format['QA'][i],BND2.format['QA'][i],BND1.format['QR'][i],BND2.format['QR'][i])
            if self.format['AB'][i] == None:
                self.format['AB'][i] = '.'
            if BND1.format['AB'][i] != '.' and BND2.format['AB'][i] != '.':
                AB1.append(float(BND1.format['AB'][i]))
                AB2.append(float(BND2.format['AB'][i]))
            i+=1
        corrcoef, pvalue = stats.pearsonr(AB1,AB2)
        self.stats["pearsonr"] = corrcoef

        #print corrcoef,pvalue


###########################################################################

# return the min_distance size of two intervals
def min_distance(BND1_start,BND1_end,BND2_start,BND2_end,overlap_space):
    "return the minium distance of two breakpoints"    
    ol_dis = min(BND2_end - BND1_start, BND1_end - BND2_start)
    #al_dis = max(BND2_end - BND1_start, BND1_end - BND2_start)
    return ol_dis


# check the orientation combination of breakpoints (0-abnormal,2-inverted,1-normal)
def check_strand(BND1_ori,BND2_ori):
    if BND1_ori == BND2_ori:
        return 0
    elif (BND1_ori == "++" and BND2_ori == "--") or (BND1_ori == "--" and BND2_ori == "++"):
        return 2
    elif (BND1_ori == "+-" and BND2_ori == "-+") or (BND1_ori == "-+" and BND2_ori == "+-"):
        return 1
    else:
        return 0


# check the order of + - reads
def check_order(oriA,oriB,posA,posB):
    "check the strand order on duplicated region"
    if posA < posB:
        if oriA == "-" and oriB == "+":
            return 1
        else:
            return 0
    else:
        if oriB == "-" and oriA == "+":
            return 1
        else:
            return 0


def check_pattern(BND1,BND2,overlap_space):
    "check if the pattern of the BND tuple match the pattern of interspersed duplication"
    
    ol1 = min_distance(BND1.start_A,BND1.end_A,BND2.start_A,BND2.end_A,overlap_space)
    ol2 = min_distance(BND1.start_B,BND1.end_B,BND2.start_B,BND2.end_B,overlap_space)
    ori = check_strand(BND1.strand_A+BND1.strand_B,BND2.strand_A+BND2.strand_B)
    if ol2 > ol1:      # always assume closer intervals should be insert position  
        ol2 += overlap_space
    else:
        ol1 += overlap_space
    # A - inserted, B - duplicated 
    if ol1 > 0 and ol2 <= 0:
        if ori  and check_order(BND1.strand_B,BND2.strand_B,BND1.start_B,BND2.start_B):
            hit = interspersed_dup(ori)
            hit.bnd1 = BND1
            hit.bnd2 = BND2
            hit.chr_ins = BND1.chrom_A
            hit.chr_dup = BND1.chrom_B
            hit.dup_start(BND1.start_B,BND2.start_B)
            hit.dup_end(BND1.end_B,BND2.end_B)
            hit.dup_insert(BND1.start_A,BND2.start_A,BND1.end_A,BND2.end_A)
            return hit	
    elif ol2 > 0 and ol1 <= 0:
        if ori and check_order(BND1.strand_A,BND2.strand_A,BND1.start_A,BND2.start_A):
            hit = interspersed_dup(ori)
            hit.bnd1 = BND1
            hit.bnd2 = BND2
            hit.chr_ins = BND1.chrom_B
            hit.chr_dup = BND1.chrom_A
            hit.dup_start(BND1.start_A,BND2.start_A)
            hit.dup_end(BND1.end_A,BND2.end_A)
            hit.dup_insert(BND1.start_B,BND2.start_B,BND1.end_B,BND2.end_B)
            return hit
    
    return None

# generate string to print in output file
def print_string(candidate):
    "Make the string to print out"
    # DUP_chrom\tDUP_start\tDUP_end\tINS_chrom\tINS_pos\tBND_1\tBND_2\tTYPE\tK_stat\tsamples\n
    s1 ='\t'.join([str(candidate.chr_dup)+':'+ \
        str(candidate.start)+'-'+ \
        str(candidate.end), \
        str(candidate.chr_ins)+':'+ \
        str(candidate.insert), \
        candidate.bnd1.id, \
        #str(candidate.sv1_NSAMP), \
        candidate.bnd2.id, \
        #str(candidate.sv2_NSAMP), \
        candidate.dup_type, \
        str(candidate.stats["pearsonr"]), \
        (str(candidate.gt_agreement[0][0])+';'+str(candidate.gt_agreement[0][1])+';'+str(candidate.gt_agreement[1][0])+';'+str(candidate.gt_agreement[1][1]))])
    s2 = '\tGT:AB\t'
    format_string = [str(x)+':'+str(y) for x,y in zip(candidate.format['GT'],candidate.format['AB']) ]
    s3 = '\t'.join(format_string)
    s = s1 + s2 + s3
    return s

def extract_intersect_BNDs(all_BNDs):
    "Extract intersect BND pairs from a BND list"
    intersect_BNDs = []
    #variants = all_BNDs.variant_list
    variants = all_BNDs.variant_list
    i = 0
    while i < len(variants):
        BND1 = variants[i]
        i += 1
        j = i 
        while j < len(variants):
            BND2 = variants[j]
            j += 1
            # Two on the same chromosome
            if BND1.chrom_A == BND1.chrom_B == BND2.chrom_A == BND2.chrom_B:
                if max(BND1.end_B,BND2.end_B) - min(BND1.start_A,BND2.start_A) < \
                (BND1.end_B - BND1.start_A) + (BND2.end_B - BND2.start_A):
                    intersect_BNDs.append((BND1,BND2))
            # On two chromosomes
            elif (BND1.chrom_A == BND2.chrom_A and BND1.chrom_B == BND2.chrom_B) :
                intersect_BNDs.append((BND1,BND2))
    return intersect_BNDs


#############################################
# primary function 
#############################################

# look for interspersed duplications
def inter_dup(bedpe_File,overlap_space,output,r_threshold):
    "Find interspersed duplications and print out"
    output = open(output,'w')
    # header
    output.write("DUP_pos\tINS_pos\tBND_1\tBND_2\tTYPE\tstat\tGT_table\tFORMAT\t")

    candidate = interspersed_dup(0)
    # Extract BND variants from bedpe file
    all_BNDs = VariantList("LUMPY_BND",bedpe_File)
    # header
    output.write('\t'.join(all_BNDs.sample_list))
    output.write('\n')

    # Extract BND pairs that are close together (for interspersed duplication cases, they must be intersect to each other)
    intersect_BNDs = extract_intersect_BNDs(all_BNDs)
    for BND_pairs in intersect_BNDs:
        BND1,BND2 = BND_pairs
        candidate = check_pattern(BND1,BND2,overlap_space)
        if candidate:
            candidate.check_GT(BND1,BND2,all_BNDs.sample_list)
            if candidate.stats['pearsonr'] > r_threshold:
                s = print_string(candidate)
                output.write(s)
                output.write("\n")

    output.close()
    #out2.close()

#############################################
# Main functions 
#############################################

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():

    usage = """%prog -b <file>
Interspersed Duplication
Author: Lei Chen
Description: find interspersed duplication.
Version:1.0
    """

    parser = OptionParser(usage)
                        
    parser.add_option("-b", "--bedpe", 
                      dest="bedpe",
                      help="bedpe file of structure variants calls")

    parser.add_option("-o","--output",
                      dest="output",
                      help="output file")
                           
    parser.add_option("-s", "--overlap_space", 
                      dest="overlap_space",
                      help="A space size in bp which indicates " + \
                           "how close two min_distance intervals should be" + \
                           "(default value = 0)",default = 0)

    parser.add_option("-r",'--r_threshold',
                      dest='r_threshold',
                      help='minium Pearson r for genotype correlation' + \
                      'of two BNDs in interspersed duplication (default=0)',default=0)

    #parser.add_option("-t", "--gt_test",
    #                  dest="gt_test",
    #                  help="Statical test for checking genotype agreement\n"+ \
    #                        "selections: kappa-gt kappa-binary pearsonr\n"+ \
    #                        "(cohen's kappa - use genotype data/Pearson R - use allele balance data)",default = "pearsonr")

    (opts, args) = parser.parse_args()

    #if opts.inFile is None or opts.configFile is None:
    if  opts.bedpe is None or opts.output is None:
        parser.print_help()
        print
    else:
        try:
            inter_dup(opts.bedpe,float(opts.overlap_space),opts.output,float(opts.r_threshold))
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
