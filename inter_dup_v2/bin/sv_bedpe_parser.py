"sv_bedpe_parser.py -- Parsing bedpe files and extract needed information"
import argparse 
import sys
import re

class VariantDetails:
	info_A = {}
	info_B = {}
	format = {}

	def __init__(self,bedpe_lines,column_list):
		"Constructor"
		self.chrom_A = None
		self.start_A = None
		self.end_A = None
		self.chrom_B = None
		self.start_B = None
		self.end_B = None
		self.id = None
		self.filter = None
		self.strand_A = None
		self.strand_B = None
		self.svtype = None
		self.column_parser(bedpe_lines,column_list)
		self.sort()
		self.info_A = {"NSAMP":0}
		self.info_B = {"NSAMP":0}
		self.format = {'GT':[],'AB':[],'QA':[],'QR':[]}
		if "INFO_A" in column_list: 
			self.info_parser(bedpe_lines[column_list.index("INFO_A")],"A")
		if "INFO_B" in column_list:
			self.info_parser(bedpe_lines[column_list.index("INFO_B")],"B")
		if "FORMAT" in column_list:
			self.format_parser(bedpe_lines[column_list.index("FORMAT")],bedpe_lines[column_list.index("FORMAT")+1:])

	def column_parser(self,bedpe_lines,column_list):
		"parse bedpe lines according to column contents"
		if "CHROM_A" in column_list:
			self.chrom_A = bedpe_lines[column_list.index("CHROM_A")]
		if "START_A" in column_list:
			self.start_A = int(bedpe_lines[column_list.index("START_A")])
		if "END_A" in column_list:
			self.end_A = int(bedpe_lines[column_list.index("END_A")])
		if "CHROM_B" in column_list:
			self.chrom_B = bedpe_lines[column_list.index("CHROM_B")]
		if "START_B" in column_list:
			self.start_B = int(bedpe_lines[column_list.index("START_B")])
		if "END_B" in column_list:
			self.end_B = int(bedpe_lines[column_list.index("END_B")])
		if "ID" in column_list:
			self.id = bedpe_lines[column_list.index("ID")]
		if "STRAND_A" in column_list:
			self.strand_A = bedpe_lines[column_list.index("STRAND_A")]
		if "STRAND_B" in column_list:
			self.strand_B = bedpe_lines[column_list.index("STRAND_B")]
		if "TYPE" in column_list:
			self.svtype = bedpe_lines[column_list.index("TYPE")]
		if "FILTER" in column_list:
			self.filter = bedpe_lines[column_list.index("FILTER")]

	def sort(self):
		"correct the order of two breakends"
		if self.chrom_A > self.chrom_B or (self.chrom_A == self.chrom_B and self.start_A > self.start_B):
			tempchr = self.chrom_A
			self.chrom_A = self.chrom_B
			self.chrom_B = tempchr
			temps1 = self.start_A
			self.start_A = self.start_B
			self.start_B = temps1
			temps2 = self.end_A
			self.end_A = self.end_B
			self.end_B = temps2
			tempori = self.strand_A
			self.strand_A = self.strand_B
			self.strand_B = tempori

	def info_parser(self,info_l,info_type):
		"Extract needed information from info columns"
		if "NSAMP" in info_l:
			temp_l = info_l[info_l.index("NSAMP"):]
			temp = temp_l.split(';')[0]
			if info_type == "A":
				self.info_A['NSAMP'] = temp.split('=')[1]
			elif info_type == "B":
				self.info_B['NSAMP'] = temp.split('=')[1]

	def  format_parser(self,formats_l,samples_l):
		"Extract needed information from formats columns and sample columns"
		GT_col = None
		AB_col = None
		QA_col = None
		QR_col = None
		format_list = formats_l.split(':')
		i = 0
		while i < len(format_list):
			if format_list[i] == 'GT':
				GT_col = i
			elif format_list[i] == 'AB':
				AB_col = i
			elif format_list[i] == 'QA':
				QA_col = i
			elif format_list[i] == 'QR':
				QR_col = i
			i += 1
		i = 0
		while i < len(samples_l):
			sample_format_list = samples_l[i].split(':')
			#print sample_format_list
			if GT_col is not None:
				self.format['GT'].append(sample_format_list[GT_col])
			if AB_col is not None:			
				self.format['AB'].append(sample_format_list[AB_col])
			if QA_col is not None:
				self.format['QA'].append(sample_format_list[QA_col])
			if QR_col is not None:
				self.format['QR'].append(sample_format_list[QR_col])
			i += 1		
		#print(len(self.format['GT']))
	
class  VariantList:
	sample_list = []
	column_list = []
	variant_list = []
	
	def __init__(self,var_pattern,bpfile):
		self.variantLister(var_pattern,bpfile)

	def check_chrom(self,var_detail):
		"Eliminate mitochondria variants, and variants in GL Contigs"
		if var_detail.chrom_A == 'M' or var_detail.chrom_B == 'M'\
		or "GL" in var_detail.chrom_A or "GL" in var_detail.chrom_A:
			return False
		else:
			return True

	def pass_filter(self,var_detail):
		"Eliminate variants didn't pass the filter"
		if var_detail.filter == 'PASS':
			return True
		else:
			return False
   

	def variantLister(self,var_pattern,bpfile):
		"Go through bedpe file and get a list of variants which match the given pattern"
		f = open(bpfile,'r')
		for l in f:	
			if "##" not in l:
				if "#" in l:
					columnline = l[l.index('#')+1:l.index('GTEX')]
					self.column_list = columnline.strip().split('\t')
					sampleline = l[l.index('GTEX'):]
					self.sample_list = sampleline.strip().split('\t')
				else:
					if re.search(var_pattern,l):
						varline = l.strip().split('\t')
						var_detail = VariantDetails(varline,self.column_list)
						#if self.check_chrom(var_detail) and self.pass_filter(var_detail):
						if self.check_chrom(var_detail):
							self.variant_list.append(var_detail)








		
		