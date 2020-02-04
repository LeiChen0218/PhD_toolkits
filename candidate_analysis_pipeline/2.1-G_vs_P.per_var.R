#####
# G_vs_P.per_var.R
# Given pheno file,variant geno file, trait, plotting geno vs pheno
# Format requirement for genotype file: three columns - SAMPLE,GT_DIS,GT_CONT (sample id, discrete genotype, continous genotype)
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2018-11-02

library(grid)
library(ggplot2)

# Collect + parse command line arguments
args <- commandArgs(T)

if(length(args) < 1){
  args <- c("--help")
}

if("--help" %in% args) {
  cat("G_vs_P.per_var.R
      Arguments:
      --geno=var_geno_file   - three columns, SAMPLE,GT_DIS,GT_CONT
      --pheno=phenotype_file   - ped format
      --trait=trait   - string
      --assoc=association_file    - association results (TRAIT column required)
      --out=output_dir - the directory to put the plots
      --help              - print this text
      ")
  
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$geno) | is.null(argsL$pheno) | is.null(argsL$trait) | is.null(argsL$assoc) | is.null(argsL$out)){
  print("missing arguements!")
  q(save = "no")
}

## test
#setwd('/Users/leichen/Desktop/Lab/general_plotting_finmetseq_sv_assoc/')

# Loading data
geno <- read.table(argsL$geno, header = T,colClasses = c("character","factor","numeric"),na.strings =".")
#geno <- read.table('candidates/CNV_chr12_18570074_18572173/data/genotype.t.list',header = T,colClasses = c("character","factor","numeric"))
pheno <- read.table(argsL$pheno, header = T, comment.char = "")
#pheno <- read.table('/Users/leichen/Desktop/Lab/finmetseq_5k_wgs_201712_201807/analysis/qt_finnseq_20161017.5k_ids.ped',header = T, comment.char = "")
trait <- argsL$trait
#trait <- "XS_VLDL_C_combined"
assoc <- read.table(argsL$assoc, header = T,colClasses = "character")
#assoc <- read.table('candidates/CNV_chr12_18570074_18572173/data/assoc.results.all_traits.sorted.txt',header = T,colClasses = "character")
output <- argsL$out
#output <- "candidates/CNV_chr12_18570074_18572173/plots"

df <- merge(geno, pheno,by.x = "SAMPLE",by.y = "IND_ID")
cand <- assoc[assoc$TRAIT==trait,]
string <- paste(list(cand),sep = "\n")

# Making plot
# function for number of observations (XXX)
give.n <- function(x){
  return(c(y = median(x)+0.3, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
# Wrap text (Richie Cotton,2010)
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

library(gtools)

out_p1 <- paste(output,"geno_vs_pheno.discrete.png",sep="/")

if("0/0" %in% df$GT_DIS){
  df$temp = df$GT_DIS
}else {
  df$temp <- factor(df$GT_DIS,levels = unique(df$GT_DIS[order(as.numeric(as.character(df$GT_DIS)))]))
}

ggplot(df,aes_string(x="temp",y=trait,fill="GT_DIS"))+
  geom_boxplot()+stat_summary(fun.data = give.n, geom = "text", fun.y = median,colour = "black")+
  scale_fill_discrete(name="Genotype/\nCopyNumber")+
  xlab("Genotype/\nCopyNumber")+
  ggtitle(wrapper(string))
ggsave(out_p1)

out_p2 <- paste(output,"geno_vs_pheno.continous.png",sep="/")
ggplot(df,aes_string(x="GT_CONT",y=trait,color="GT_DIS"))+
  geom_point()+xlab("Allele Balance/\nFractional Copy Number")+
  ggtitle(wrapper(string))
ggsave(out_p2)