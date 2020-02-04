#####
# plot_cn_profile.per_var.R
# Given cn file, signal CNV info, signal CNV genotype, plotting copy number profile
# Format requirement for genotype file: three columns - SAMPLE,GT_DIS,GT_CONT (sample id, discrete genotype, continous genotype)
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2018-12-12

source('/Users/leichen/Desktop/Lab/Finmetseq_paper/3-Candidate_analysis/scripts/bin/6.1-cn.heatmap.fun.R')

# Collect + parse command line arguments
args <- commandArgs(T)

if(length(args) < 1){
  args <- c("--help")
}

if("--help" %in% args) {
  cat("plot_cn_profile.per_var.R
      Arguments:
      --geno=var_geno_file   - three columns, SAMPLE,GT_DIS,GT_CONT
      --cn=cn_file   - copy number matrix (format similar to:http://software.broadinstitute.org/cancer/software/genepattern/tutorial/linkedFiles/mynah.sorted.cn) 
      --sig=signal   - string
      --info=var_info_file    - SV signal info (col:VAR,CHR,POS,END,TYPE,LEN,AF)
      --out=output - the output directory (and prefix) of the plots
      --left_window=left_plot_boundary - (optional) the number of regions to be included in the plot left to SV [maximum size allowed in data]
      --right_window=right_plot_boundary - (optional) the number of regions to be included in the plot right to SV [maximum size allowed in data]
      --cn_sat=saturation_cn - (optional) upper limit of color saturation for copy number value [4]
      --help              - print this text
      ")
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$geno) | is.null(argsL$cn) | is.null(argsL$sig) | is.null(argsL$info) | is.null(argsL$out)){
  print("missing arguements!")
  q(save = "no")
}

############
# read in tables and call functions 

# test the function
#setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/3-Candidate_analysis/data/cand_plots/chr2_1343601_1344800/data')

cn_matrix <- read.table(argsL$cn, header = T)

#cn_matrix <-  read.table('rep100sample.cn',header = T)
#cn_heatmap_fun(cn_matrix,start=start, end=end,out=out)
#print(argsL$geno)
dosage <- read.table(argsL$geno, header = T)
#dosage <- read.table('genotype.t.list',header = T)
dosage_order <- get_dosage_order(cn_matrix,dosage)
#info <- read.table('../../cands/cand.var.info')
info <- read.table(argsL$info)
colnames(info) <- c("VAR","CHR","POS","END","TYPE","LEN","AF")
sig <- info[info$VAR==argsL$sig,]
#sig <- info[info$VAR=="chr2_1343601_1344800",]

### get the start and end by given left window numbers and right window numbers

t=unlist(strsplit(as.character(cn_matrix$SNP[1]), "[:-]", perl = T))
window_size <- as.numeric(t[3]) - as.numeric(t[2])

# select the region to plot
if(is.null(argsL$left_window)){
  start = min(cn_matrix$PhysicalPosition)
}else{
  left_window_size = as.numeric(argsL$left_window)*window_size
  start = sig$POS - left_window_size
  }
if(is.null(argsL$right_window)){
  end = max(cn_matrix$PhysicalPosition)
}else{
  right_window_size = as.numeric(argsL$right_window)*window_size
  end = sig$END + right_window_size
}

#out="/Users/leichen/Desktop/Lab/general_plotting_finmetseq_sv_assoc/candidates/44368/plots/w100bp"
out=paste(argsL$out,"cn",sig$VAR,"rep100samp_order_by_GT",paste(sig$CHR,start,end,sep="_"),sep=".")
cn_heatmap_fun(cn_matrix,order_ref=dosage_order,start=start, end=end,cn_sat=argsL$cn_sat,sig=sig,out=out)
out2=paste(argsL$out,"cn",sig$VAR,"rep100samp_order_by_clustering",paste(sig$CHR,start,end,sep="_"),sep=".")
cn_heatmap_fun(cn_matrix=cn_matrix,start=start, end=end,cn_sat=argsL$cn_sat,sig=sig,out=out2)
