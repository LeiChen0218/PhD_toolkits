#! /usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
import numpy as np
import pandas as pd
from scipy import stats
from collections import Counter, defaultdict, namedtuple
import statsmodels.formula.api as smf
from operator import itemgetter
import warnings
from sv_vcf_parser import Vcf, Variant, Genotype
import os

CN_rec = namedtuple ('CN_rec', 'var_id sample svtype svlen AF GT CN AB log_len log2r')

class idup_variant:
    
    def __init__(self, var_line,samplelist):
        self.dup_pos = var_line[0]
        self.chrom = var_line[0].split(':')[0]
        self.ins_pos = var_line[1]
        self.BND1 = var_line[2]
        self.BND2 = var_line[3]
        self.ins_ori = var_line[4]
        self.format = {}
        self.sample_list = samplelist
        self.active_formats = []
        self.parse_format(var_line[8:])
        self.freq_type = None

    def parse_format(self,format_line):
        i = 0
        while i < len(self.sample_list):
            self.format[self.sample_list[i]] = {'GT':format_line[i].split(':')[0],'AB':format_line[i].split(':')[1],'CN':None}
            i += 1
        self.active_formats.append('GT')
        self.active_formats.append('AB')

    def add_cn(self, cn_dic):
        i = 0
        while i < len(self.sample_list):
            if self.sample_list[i] in cn_dic[self.dup_pos].keys():
                self.format[self.sample_list[i]]['CN'] = float(cn_dic[self.dup_pos][self.sample_list[i]])*2
            else:
                self.format[self.sample_list[i]]['CN'] = '.'
            i += 1
        self.active_formats.append('CN')

# http://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def to_bnd_strings(var, fixed_gts):

    old_type = var.info['SVTYPE']
    old_id = var.var_id
    old_pos = var.pos
    old_end = var.info['END']
    old_ciend = var.info['CIEND']
    old_cipos = var.info['CIPOS']
    old_cipos95 = var.info['CIPOS95']
    old_ciend95 = var.info['CIEND95']


    #for both ends
    var.info['SVTYPE'] = 'BND'
    var.info['EVENT'] = old_id
    del var.info['SVLEN']
    del var.info['END']

    #var1
    var.var_id = old_id + "_1"
    var.info['MATEID'] = old_id + "_2"
    if old_type == 'DEL':
        var.alt = 'N[%s:%s[' % (var.chrom, old_end)
    else:
        var.alt = ']%s:%s]N' % (var.chrom, old_end)
    var1=var.get_var_string(fixed_gts)

    #var2
    var.var_id = old_id + "_2"
    var.info['MATEID'] = old_id + "_1"
    var.info['CIPOS'] = old_ciend
    var.info['CIEND'] = old_cipos
    var.info['CIPOS95'] = old_ciend95
    var.info['CIEND95'] = old_cipos95
    var.pos = old_end
    var.info['SECONDARY'] = True
    if old_type == 'DEL':
        var.alt = ']%s:%s]N' % (var.chrom, old_pos)
    else:
        var.alt = 'N[%s:%s[' % (var.chrom, old_pos)
    var2=var.get_var_string(fixed_gts)
    return var1, var2


def lowQuantile(xx):
    return np.percentile(xx,2.5)

def highQuantile(xx):
    return np.percentile(xx,97.5)

def lld(xx, mean, sd):
    ll = 1 / sd * math.exp(-(xx-mean) * (xx-mean) / (2*sd*sd))
    return ll

def calc_params(vcf_path):

    tSet = list()
    epsilon=0.1
    header=[]
    
    in_header = True
    vcf = Vcf()
    if vcf_path.endswith('.gz'):
        vcf_file = gzip.open(vcf_path, 'rb')
    else:
        vcf_file = open(vcf_path, 'r')

    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line)
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                    in_header = False
                    vcf.add_header(header)
                continue
        else:
            v = line.rstrip().split('\t')
            info = v[7].split(';')
            svtype = None
            for x in info:
                if x.startswith('SVTYPE='):
                    svtype = x.split('=')[1]
                    break

            if svtype not in ['DEL', 'DUP'] or v[0]=="X" or v[0]=="Y":
                continue

            var = Variant(v, vcf)
    
            for sample in vcf_samples:
                sample_genotype = var.genotype(sample)
                if sample_genotype.get_format('GT') != './.':
                    log2r = math.log((float(sample_genotype.get_format('CN'))+ epsilon)/2,2)  #to avoid log(0)
                    tSet.append(CN_rec(var.var_id, sample, var.info['SVTYPE'], abs(float(var.info['SVLEN'])), var.info['AF'],
                        sample_genotype.get_format('GT'),  sample_genotype.get_format('CN'), sample_genotype.get_format('AB'), math.log(abs(float(var.info['SVLEN']))), log2r))

    df=pd.DataFrame(tSet, columns=CN_rec._fields)
    #exclude from training data, DELs and DUPs with CN in the tails of the distribution
    df['q_low']=df.groupby(['sample', 'svtype', 'GT'])['log2r'].transform(lowQuantile)
    df['q_high']=df.groupby(['sample', 'svtype', 'GT'])['log2r'].transform(highQuantile)
    df=df[(df.log2r>=df.q_low) & (df.log2r<=df.q_high)]
    #df.to_csv('./train.csv')

    #adjust copy number for small deletions (<1kb), no strong relationship b/w cn and size for dups evident so far
    small_het_dels = df[(df.svtype=="DEL") & (df.GT=="0/1") & (df.svlen<1000) & (df.svlen>=50)]
    small_hom_dels = df[(df.svtype=="DEL") & (df.GT=="1/1") & (df.svlen<1000) & (df.svlen>=50)]
    het_del_mean=np.mean(df[(df.svlen>1000) & (df.GT=="0/1") & (df.svtype=="DEL")]['log2r'])
    hom_del_mean=np.mean(df[(df.svlen>1000) & (df.GT=="1/1") & (df.svtype=="DEL")]['log2r'])
    small_het_dels['offset']=small_het_dels['log2r']-het_del_mean
    small_hom_dels['offset']=small_hom_dels['log2r']-hom_del_mean
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        hom_del_fit=smf.ols('offset~log_len',small_hom_dels).fit()
        het_del_fit=smf.ols('offset~log_len',small_het_dels).fit()
        #print hom_del_fit.summary()
        #print het_del_fit.summary()
        small_hom_dels['log2r_adj'] = small_hom_dels['log2r'] - hom_del_fit.predict(small_hom_dels)
        small_het_dels['log2r_adj'] = small_het_dels['log2r'] - het_del_fit.predict(small_het_dels)

    small_dels=small_hom_dels.append(small_het_dels)
    small_dels=small_dels[['var_id', 'sample', 'svtype', 'svlen', 'AF', 'GT', 'CN', 'log_len', 'log2r', 'q_low', 'q_high', 'log2r_adj']]

    # dels of length<100 bp are excluded here
    df1=df[(df.svtype!="DEL") | (df.GT=="0/0") | (df.svlen>=1000)]
    df1['log2r_adj']=df1['log2r']
    df1=df1.append(small_dels)

    params=df1.groupby(['sample', 'svtype', 'GT'])['log2r_adj'].aggregate([np.mean,np.var, len]).reset_index()
    params=pd.pivot_table(params, index=['sample', 'svtype'], columns='GT', values=['mean', 'var', 'len']).reset_index()    
    params.columns=['sample', 'svtype', 'mean0', 'mean1', 'mean2', 'var0', 'var1', 'var2', 'len0', 'len1', 'len2']
    params['std_pooled']=np.sqrt((params['var0']*params['len0']+params['var1']*params['len1']+params['var2']*params['len2'])/(params['len0']+params['len1']+params['len2']))
    #params.to_csv('./params.csv')
    return (params, het_del_fit, hom_del_fit)


def rd_support_nb(temp, p_cnv):
    tr = pd.DataFrame({'p0' : [1.0, 0.1, 0.0], 'p1' : [0.0, 0.7, 0.25], 'p2' : [0.0, 0.2, 0.75], 'GT' : ["0/0", "0/1", "1/1"]})
    temp = pd.merge(temp, tr, on='GT', how='left')
    temp['p_mix'] = temp['lld0'] * temp['p0'] + temp['lld1'] * temp['p1'] + temp['lld2'] * temp['p2']
    return np.log(p_cnv)+np.sum(np.log(temp['p_mix'])) > np.log(1-p_cnv)+np.sum(np.log(temp['lld0']))
   

def has_rd_support_by_nb(test_set, params, p_cnv = 0.5):
    svtype=test_set['svtype'][0]
    svlen=test_set['svlen'][0]
    log_len=test_set['log_len'][0]
    
    params1=params.copy()
    params1['mean1_adj'] = params1['mean1']
    params1['mean2_adj'] = params1['mean2']

    v0=test_set[test_set.GT=="0/0"]['log2r'].values
    v1=test_set[test_set.GT=="0/1"]['log2r'].values
    v2=test_set[test_set.GT=="1/1"]['log2r'].values

    if len(v0)>0:
        med0=np.median(v0)
    else:
        if len(v1)>0:
            med0=med1=np.median(v1)
        elif len(v2)>0:
            med0=med1=med2=np.median(v2)
        else:
            return False

    if len(v1)>0:
        med1=np.median(v1)
    else:
        med1=med0
    if len(v2)>0:
        med2=np.median(v2)
    else:
        med2=med1

    if svtype=='DEL' and ( med1>med0 or med2>med0 ):
        return False
    elif svtype=='DUP' and (med1<med0 or med2<med0):
        return False

    mm=pd.merge(test_set, params1, how='left')

    mm['lld0'] = mm.apply(lambda row:lld(row["log2r"], row["mean0"],row["std_pooled"]), axis=1)
    mm['lld1'] = mm.apply(lambda row:lld(row["log2r"], row["mean1_adj"],row["std_pooled"]), axis=1)
    mm['lld2'] = mm.apply(lambda row:lld(row["log2r"], row["mean2_adj"],row["std_pooled"]), axis=1)
   
    return  rd_support_nb(mm, p_cnv)

def get_sv_len(dup_pos):

	start = float(dup_pos.split(':')[1].split('-')[0])
	end = float(dup_pos.split(':')[1].split('-')[1])

	return abs(end - start)

def get_af(var_format):

	AA_count = 0
	RA_count = 0
	GT_list = []

	for s in var_format.keys():
		GT_list.append(var_format[s]['GT'])

	for gt in GT_list:
		if gt == '1/1':
			AA_count += 2.0
		elif gt == '0/1':
			AA_count += 1.0
			RA_count += 1.0
		elif gt == '0/0':
			RA_count += 2.0

	return AA_count/(RA_count+AA_count)

def load_df(var, exclude, sex):
    
    epsilon=0.1
    test_set = list()

    for s in var.sample_list:
        if s in exclude:
            continue
        cn = var.format[s]['CN']
        if (var.chrom == 'X' or var.chrom == 'Y') and sex[s] == 1:
            cn=str(float(cn)*2)
        log2r = math.log((float(cn)+epsilon)/2, 2)  # to avoid log(0)
        test_set.append(CN_rec(var.dup_pos, s, 'DUP', get_sv_len(var.dup_pos), get_af(var.format),
            var.format[s]['GT'],  cn , var.format[s]['AB'], math.log(get_sv_len(var.dup_pos)), log2r))

    test_set = pd.DataFrame(data = test_set, columns=CN_rec._fields)
    return test_set

# test for read depth support of low frequency variants
def has_low_freq_depth_support(test_set, mad_threshold=2, absolute_cn_diff=0.5):

    mad_quorum = 0.5 # this fraction of the pos. genotyped results must meet the mad_threshold
    
    hom_ref_cn=test_set[test_set.GT=="0/0"]['CN'].values.astype(float)
    hom_het_alt_cn=test_set[(test_set.GT=="0/1") | (test_set.GT=="1/1")]['CN'].values.astype(float)

    if len(hom_ref_cn) > 0:
        cn_median = np.median(hom_ref_cn)
        cn_mad = mad(hom_ref_cn)
    else:
        cn_median = None
        cn_mad = None

    # bail after writing out diagnostic info, if no ref samples or all ref samples
    if (len(hom_ref_cn) == 0 or len(hom_het_alt_cn) == 0):
        return False

    # tally up the pos. genotyped samples meeting the mad_threshold

    resid=hom_het_alt_cn-cn_median
    if test_set['svtype'][0]=='DEL':
        resid=-resid
    
    resid=resid[(resid > (cn_mad * mad_threshold) ) & (resid>absolute_cn_diff)]

    if float(len(resid))/len(hom_het_alt_cn)>mad_quorum:
        return True
    else:
        return False

# test whether variant has read depth support by regression
def has_high_freq_depth_support(df, slope_threshold, rsquared_threshold):
    
    rd = df[[ 'AB', 'CN']][df['AB']!='.'].values.astype(float)
    if len(np.unique(rd[0,:])) > 1 and len(np.unique(rd[1,:])) > 1:
        
        (slope, intercept, r_value, p_value, std_err) = stats.linregress(rd)

        if (slope < slope_threshold or r_value*r_value < rsquared_threshold):
            return False
        return True
    return False

def has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps, mad_threshold=2, absolute_cn_diff=0.5):

    min_pos_samps_for_regression=10
    if num_pos_samps>min_pos_samps_for_regression:
        return has_high_freq_depth_support(df, slope_threshold, rsquared_threshold)
    else:
        return has_low_freq_depth_support(df, mad_threshold, absolute_cn_diff)
    return False

def has_rd_support_hybrid(df, params, p_cnv, slope_threshold, rsquared_threshold, num_pos_samps):
  
    hybrid_support=False
    nb_support=has_rd_support_by_nb(df, params, p_cnv)
    ls_support=has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps)

    if nb_support and ls_support:
        hybrid_support=True
    elif nb_support and has_rd_support_by_ls(df, 2*slope_threshold, 2*rsquared_threshold, num_pos_samps, 2, 0.75):
        hybrid_support=True
    elif ls_support and has_rd_support_by_nb(df, params, 0.2*p_cnv):
        hybrid_support=True
    return [ls_support, nb_support, hybrid_support]


# primary function
def sv_classify(cn_dic,
                idup_var,
                gender_file, 
                exclude_file, 
                slope_threshold, 
                rsquared_threshold, 
                p_cnv, 
                params,
                method):

    #out_temp = open('valid.txt','w')
    sample_column_index = 0
    sample_list = []
    gender = {}
    num_pos_samps = 0
    # read sample genders
    for line in gender_file:
        v = line.rstrip().split('\t')
        sample_short = '-'.join(v[0].split('-')[0:2])
        gender[sample_short] = int(v[1])

    exclude = []
    if exclude_file is not None:
        for line in exclude_file:
            exclude.append(line.rstrip())
    exclude = ['GTEX-NPJ8','GTEX-S7PM']

    valid_BNDs = {}
    for line in idup_var:
        if 'GTEX' in line:
            sample_column_index = line.index('GTEX')
            sample_list = line[sample_column_index:].strip().split()
        else:
            var = idup_variant(line.strip().split('\t'),sample_list)
            var.add_cn(cn_dic)
            # count the number of positively genotyped samples
            num_pos_samps = 0;
            for s in var.sample_list:
                if s in exclude:
                    continue
                if var.format[s]['GT'] not in ["./.", "0/0"]:
                    num_pos_samps += 1

        nb_support = False
        ls_support = False
        hybrid_support = False
        has_rd_support = False

        if num_pos_samps == 0:
            continue
        else:
            df=load_df(var, exclude, gender)
            if method=='large_sample':
                ls_support = has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps)
                has_rd_support=ls_support
                
            elif method=='naive_bayes':
                nb_support = has_rd_support_by_nb(df, params, p_cnv)
                has_rd_support=nb_support
            elif method=='hybrid':
                ls_support, nb_support, hybrid_support = has_rd_support_hybrid(df, params, p_cnv, slope_threshold, rsquared_threshold, num_pos_samps)
                has_rd_support=hybrid_support

            if has_rd_support:
                #out_temp.write(str(var.dup_pos)+'\n')
                valid_BNDs[var.BND1] = var
                valid_BNDs[var.BND2] = var
                
    gender_file.close()
    idup_var.close()
    if exclude_file is not None:
        exclude_file.close()

    return valid_BNDs

def parse_cn_file(cn_file):
    "Parse the results file from CNVnator"
    cn_dic = {}
    dirs = os.listdir(cn_file)
    for filename in dirs:
        if 'txt' in filename:
            abpath_file = os.path.join(cn_file,filename)
            c_f = open(abpath_file,'r')
        else:
            continue
        for l in c_f:
            line = l.strip().split()
            if len(line) < 4:
                continue
            if line[1] not in cn_dic.keys():
                cn_dic[line[1]] = {}
            
            sample = filename.split('.')[0]
            cn_dic[line[1]][sample] = line[3]
    
    return cn_dic

def set_BNDTYPE(var,idup):

    var.set_info('BNDTYPE', 'iDUP_'+idup.ins_ori)

    return var

def set_RELATEDBND(bnd_id,var,idup):
    mated = var.get_info('MATEID')
    if bnd_id == idup.BND1:
        var.set_info('MATEID',','.join((mated,idup.BND2+'_1',idup.BND2+'_2')))
    elif bnd_id == idup.BND2:
        var.set_info('MATEID',','.join((mated,idup.BND1+'_1',idup.BND1+'_2')))
    return var

def modify_info_in_vcf(vcf_in,vcf_out,valid_BNDs):
    in_header = True
    header = []
    vcf = Vcf()
    BNDTYPE_description = "Variant type after BND reclassification"
    DUPPOS_description = "Duplication region position"
    RELATEDBND_description = "Related BND(s) in reclassification"
    INSPOS_description = "Insertion position"
    # read input VCF
    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line) 
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_info('BNDTYPE', 1, 'String', BNDTYPE_description)
                vcf.add_info('DUPPOS', 1, 'String', DUPPOS_description)
                vcf.add_info('INSPOS', 1, 'String', INSPOS_description)
                #vcf.add_info('RELATEDBND', '.', 'String', RELATEDBND_description)
                vcf_out.write(vcf.get_header(include_samples=True)+'\n')
        else:
            if 'LUMPY_BND' in line:
                var_list = line.split('\t')
                bnd_id = re.findall(r'LUMPY_BND_[0-9]+',var_list[2])[0]

                if bnd_id in valid_BNDs.keys():
                    sys.stderr.write("Matched original BND\n")
                    var = Variant(var_list, vcf)
                    idup = valid_BNDs[bnd_id]
                    # Add info
                    for existed_info in var_list[7].split(';'):
                        if '=' in existed_info:
                            e_info_id = existed_info.split('=')[0]
                            e_info_value = existed_info.split('=')[1]
                            var.set_info(e_info_id,e_info_value)
                    var = set_BNDTYPE(var,idup)
                    var = set_RELATEDBND(bnd_id,var,idup)
                    var.set_info('DUPPOS',idup.dup_pos)
                    var.set_info('INSPOS',idup.ins_pos)
                    # Modify GT,AB and add CN
                    for s in var.sample_list:
                        var.gts[s].set_format('GT',idup.format[s]['GT'])
                        var.gts[s].set_format('CN',idup.format[s]['CN'])
                        if idup.format[s]['AB'] != '.':
                            var.gts[s].set_format('AB',float(idup.format[s]['AB']))
                    vcf_out.write(var.get_var_string())
                else:
                    vcf_out.write(line)
            else:
                vcf_out.write(line)

    vcf_out.close()
    vcf_in.close()


def run_reclassifier(vcf_file,
                     idup_var,
                     cn_file,
                     vcf_out,
                     sex_file,
                     exclude_list,
                     slope_threshold,
                     rsquared_threshold,
                     training_data,
                     method):

    ae_dict = None
    params = None
    het_del_fit = None
    hom_del_fit = None
    cn_dic = None
    valid_BNDs = {}
    p_cnv=0.5       # prior probability that CNV is real
    
    if(method!="large_sample"):
          sys.stderr.write("calculating parameters\n")
          #calculate per-sample CN profiles on training set
          [params, het_del_fit, hom_del_fit]=calc_params(training_data)

    sys.stderr.write("copy number file loading\n")

    cn_dic = parse_cn_file(cn_file)
    valid_BNDs = sv_classify(cn_dic,
                             idup_var,
                             sex_file,
                             exclude_list,
                             slope_threshold,
                             rsquared_threshold,
                             p_cnv,
                             params,
                             method)
    sys.stderr.write("copy number test has been done\nvalid variants:"+str(len(valid_BNDs.keys())/2)+'\n')
    modify_info_in_vcf(vcf_file,vcf_out,valid_BNDs)


def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-v', '--variants',metavar='<FILE>',dest='idup_var',type=argparse.FileType('r'),default=None,required=True,help='interspersed duplication variants file')
    parser.add_argument('-c', '--copynumber',metavar='<STRING>',dest='cn_file',type=str,default=None,required=True,help='Path to CNVnator results file')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-g', '--gender', metavar='<FILE>', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.add_argument('-e', '--exclude', metavar='<FILE>', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.add_argument('-s', '--slope_threshold', metavar='<FLOAT>', dest='slope_threshold', type=float, default=1.0, help='minimum slope absolute value of regression line to classify as DEL or DUP[1.0]')
    parser.add_argument('-r', '--rsquared_threshold', metavar='<FLOAT>', dest='rsquared_threshold', type=float, default=0.2, help='minimum R^2 correlation value of regression line to classify as DEL or DUP [0.2], for large sample reclassification')
    parser.add_argument('-t', '--tSet', metavar='<STRING>', dest='tSet', type=str, default=None, required=False, help='high quality deletions & duplications training dataset[vcf], required by naive Bayes reclassification')
    parser.add_argument('-m', '--method', metavar='<STRING>', dest='method', type=str, default="large_sample", required=False, help='reclassification method, one of (large_sample, naive_bayes, hybrid)', choices=['large_sample', 'naive_bayes', 'hybrid'])
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'reclassify interspersed duplications based on read depth information'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    if args.tSet is None:
        if args.method!="large_sample":
            sys.stderr.write("Training data required for naive Bayes or hybrid classifiers\n")
            parser.print_help()
            sys.exit(1)
    run_reclassifier(args.vcf_in, 
    	             args.idup_var,
    	             args.cn_file, 
    	             args.vcf_out, 
    	             args.gender,
    	             args.exclude, 
    	             args.slope_threshold, 
    	             args.rsquared_threshold, 
    	             args.tSet, 
    	             args.method)

if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
