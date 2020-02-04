args <- commandArgs(trailingOnly = T)
if(length(args) < 1)
  stop("input file missing")
epact_results <- args[1]

if( length(grep(".gz",epact_results)) > 0 ){
  epact<-read.table(gzfile(epact_results),header = T,comment.char = "")
}     else{
  epact<-read.table(epact_results,header = T,comment.char = "")}

#setwd('/Users/leichen/Desktop/Lab/finmetseq_5k_wgs_201712_201807/outputs/S_krea_combined')
#epact <- read.table('finmetseq.5k.S_krea_combined.emmax.kin.lm.assoc.epacts.gz', header = T,comment.char = "")
library(grid)
library(ggplot2)
library(qqman)

make_manhattan_plot <- function(df){
  selected_col <- c("MARKER_ID","X.CHROM","BEG","PVALUE","BETA")
  gwas_format <- df[selected_col]
  colnames(gwas_format) <- c("SNP","CHR","BP","P","zscore")
  if("X" %in% gwas_format$CHR)
     gwas_format$CHR <- gsub("X","23",gwas_format$CHR)
  if("Y" %in% gwas_format$CHR)
     gwas_format$CHR <- gsub("Y","24",gwas_format$CHR)
  gwas_format$CHR <- as.numeric(gwas_format$CHR)
  gwas_format_subset <- na.omit(gwas_format)
  manhattan(gwas_format_subset)
  # Add p_value thredhold line
  abline(h=3.38e-8,col = 'gray')
  return(0)
}

make_QQ_plot <- function(df){
  logp <- -log10(na.omit(df$PVALUE))
  exp_p <- seq(from=1/length(logp), to= 1, by = 1/length(logp))
  qq.df <- as.data.frame(cbind(sort(logp, decreasing = T), -log10(exp_p)))
  plot(V1~V2, data=qq.df, xlab="Expacted ordered -log10 Pvalue", ylab = "ordered -log10 Pvalue")
  abline(a=0,b=1,col="red")
}
# 
# # 

if( length(grep(".gz",epact_results)) > 0 ){
  output <- gsub("gz","mh.pdf",epact_results)
  output1 <- gsub("gz","qq.pdf",epact_results)
  output2 <- gsub("gz","qq.autosome.pdf",epact_results)
}     else{
  output <- paste(epact_results,"mh.pdf", sep=".")
  output1 <- paste(epact_results,"qq.pdf", sep=".")
  output2 <- paste(epact_results,"qq.autosome.pdf", sep=".")
}


pdf(output,width = 12,height = 8)
 make_manhattan_plot(epact)
dev.off()
# # 
 
pdf(output1,width = 8,height = 8)
make_QQ_plot(epact)
 dev.off()

pdf(output2,width = 8,height = 8)
make_QQ_plot(epact[epact$X.CHROM != "X" & epact$X.CHROM !="Y",])
dev.off()

