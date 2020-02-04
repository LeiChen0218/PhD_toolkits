#####
# Phewide_Mplot.per_var.R
# Given pheno group file,association results file, plotting trait-wise p-value plot 
# Format requirement for association file: TRAIT	VAR	CHR	POS	PVALUE	BETA
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2018-11-05

library(grid)
library(ggplot2)

# Collect + parse command line arguments
args <- commandArgs(T)

if(length(args) < 1){
  args <- c("--help")
}

if("--help" %in% args) {
  cat("Phewide_Mplot.per_var.R
      Arguments:
      --group=phenotype_group_file  - csv format, no header
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

if(is.null(argsL$group) | is.null(argsL$assoc) | is.null(argsL$out)){
  print("missing arguements!")
  q(save = "no")
}

# test script
#setwd('/Users/leichen/Desktop/Lab/general_plotting_finmetseq_sv_assoc/')

pheno.group <- read.csv(argsL$group, header = F)
#pheno.group <- read.csv("/Users/leichen/Desktop/Lab/coldspot/gtex_eqtl/portal_499s_44t/counts/large_ss_combined_pheno_201610.csv",header = F)
colnames(pheno.group) <- c("TRAIT","group")

assoc.df = read.table(argsL$assoc,header = T)
#assoc.df = assoc <- read.table('candidates/CNV_chr16_69982813_69983880/data/assoc.results.all_traits.sorted.txt',header = T)
assoc.df$logP <- -log10(assoc.df$PVALUE)
assoc.df$trait <- gsub('_rn','',assoc.df$TRAIT) # modify the trait name

df <- merge(assoc.df,pheno.group, by.x="trait" ,by.y="TRAIT")
assoc.sorted.df <- df[order(df$PVALUE),]
assoc.sorted.df$TRAIT <- factor(assoc.sorted.df$TRAIT,levels = assoc.sorted.df$TRAIT[order(assoc.sorted.df$PVALUE)])

library("ggthemes")
out <- paste(argsL$out,"all_traits_p.png",sep="/")
ggplot(assoc.sorted.df,aes(x=TRAIT, y=logP, fill=group)) + 
  geom_bar(stat="identity") +
  scale_fill_gdocs()+ylab('-log10 p-value')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 4))+
  ggtitle(apply(assoc.sorted.df[1,c(1:2,5)],1,paste,collapse="\n"))
ggsave(out,width = 12,height = 8)

