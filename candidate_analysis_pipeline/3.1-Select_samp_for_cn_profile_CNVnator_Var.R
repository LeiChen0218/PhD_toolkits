#####
# Select_samp_for_cn_profile_CNVnator.R
# Given genotype file, select samples to use in the later copy number profile plot
# Format requirement : column1:sample id column2:copy number
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2019-03-08
# Collect + parse command line arguments
args <- commandArgs(T)

if(length(args) < 1){
  args <- c("--help")
}

if("--help" %in% args) {
  cat("Phewide_Mplot.per_var.R
      Arguments:
      --geno=genotype_file  - column1:sample_id column2:genotype
      --stat=caller_stat - merge6 output format
      --out=output_dir - the directory to put the sample list
      --help              - print this text
      ")
  
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$geno) | is.null(argsL$out)){
  print("missing arguements!")
  q(save = "no")
}

# test script
#setwd('/Users/leichen/Desktop/Lab/general_plotting_finmetseq_sv_assoc/')
#geno <- read.table('candidates/103068//data/genotype.discrete.t.list',header = T, colClasses = c("character","factor"))
geno <- read.table(argsL$geno,header = T, colClasses = c("character","factor") )
stat <- read.table(argsL$stat, header = T)
colnames(geno) <- c("ID","GT")
var_stat <- stat[1,]

####
# 

# if variant is rare, just take top 20 + bottom 20 + random 60
if(var_stat$is_rare == 1){
  
  
}else{ # otherwise classify the samples into cn groups and then use the same strategy for genomestrip CNVs
  
}
  



ngroup <- var_stat$nocl






samples <- vector('character')
for(i in 1:length(geno_count)){
  sample_pool <- geno$ID[geno$GT == names(geno_count[i])]
  
  if(geno_count[i] <= ngroup)
    samples <- c(samples,sample_pool)
  else if(i == length(geno_count))
    samples <- c(samples,sample(sample_pool,100-length(samples)))
  else
    samples <- c(samples,sample(sample_pool,ngroup))
}

#double_check <- geno[geno$ID %in% samples,]
#table(double_check$GT)

#out="candidates/CNV_chr16_69982813_69983880/data"
output <- paste(argsL$out,"rep100sample.list",sep="/")
write(samples,file=output,sep="\n")
