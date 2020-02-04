#####
# ld.snp.plot.fun.R
# Given genotype corrlation file ,sv signal file, snp association results, trait, plotting ld & association signal info
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2018-12-05

library(grid)
library(ggplot2)

# major plotting function
plot_ld_w_snps <- function(sig_table,
                           trait,
                           ld,
                           assoc,
                           out_dir){
  
  
  # join the snp tables
  assoc$logP <- -log10(assoc$PVALUE)
  assoc$snp1 = paste(assoc$X.CHROM,assoc$BEG,sep = ":")
  snps <- merge(assoc,ld,by = "snp1")
  
  best_tag <- snps[snps$rsquare == max(snps$rsquare), ]
  snp_out <- paste(out_dir,"_best_tagged_snps.txt",sep="")
  write.table(best_tag,file = snp_out,quote = F,sep="\t",row.names = F)
  
  # get the row of signal sv and reformat
  sig <- sig_table[sig_table$TRAIT == trait,]
  colnames(sig) <- c("TRAIT","VAR","CHR","BEG","PVALUE","BETA")
  sig$logP <- -log10(sig$PVALUE)
  
  # plot
  neighborhood<-ggplot(snps, aes(x = BEG, y = logP,color=rsquare)) + geom_point()
  ld_plot <- neighborhood + geom_point(data = sig,color="red",shape=9,size=6) +  scale_colour_gradient(low = "green", high = "red",limits=c(0, 1))+
    annotate("text", label =sig$VAR, x = sig$BEG, y = sig$logP-0.5, size = 4, colour = "black")+
    annotate("text", label = paste("p-value = ",sig$PVALUE), x = sig$BEG, y = sig$logP-1, size = 4, colour = "black")+
    ylab("-log10(P.value)")+xlab(paste("position on chromosome",sig$CHR))+
    ggtitle(paste(trait,sig$VAR))
  
  label_best_tag_snp <- paste("max r2",best_tag[1,]$rsquare, sep="=")
  label_best_tag_snp <- paste(label_best_tag_snp,best_tag[1,]$snp1, sep = "\n")
  p <- ld_plot + annotate("text", label = label_best_tag_snp, x = best_tag[1,]$BEG, y =  best_tag[1,]$logP+0.5, size = 3, colour = "black")
  pp <- p + geom_segment(x= best_tag[1,]$BEG, xend =best_tag[1,]$BEG, y= best_tag[1,]$logP+0.3, yend = best_tag[1,]$logP, size = 0.01, arrow =  arrow(length = unit(0.01, "npc")), colour = "black")
  
  # save image file
  plot_out <- paste(out_dir,"ld_plot_",trait,"_",sig$VAR,".png",sep = "")
  print(pp)
  ggsave(plot_out,  width = 10,height = 6)
  
  # plot snps with R2 > 0 only
  if(max(snps$rsquare) <= 0.001)
    return()
  
  neighborhood<-ggplot(snps[snps$rsquare > 0.001,], aes(x = BEG, y = logP,color=rsquare)) + geom_point()
  ld_plot <- neighborhood + geom_point(data = sig,color="red",shape=9,size=6) +  scale_colour_gradient(low = "green", high = "red",limits=c(0, 1))+
    annotate("text", label =sig$VAR, x = sig$BEG, y = sig$logP-0.5, size = 4, colour = "black")+
    annotate("text", label = paste("p-value = ",sig$PVALUE), x = sig$BEG, y = sig$logP-1, size = 4, colour = "black")+
    ylab("-log10(P.value)")+xlab(paste("position on chromosome",sig$CHR))+
    ggtitle(paste(trait,sig$VAR))
  
  label_best_tag_snp <- paste("max r2",best_tag[1,]$rsquare, sep="=")
  label_best_tag_snp <- paste(label_best_tag_snp,best_tag[1,]$snp1, sep = "\n")
  p <- ld_plot + annotate("text", label = label_best_tag_snp, x = best_tag[1,]$BEG, y =  best_tag[1,]$logP+0.5, size = 3, colour = "black")
  pp <- p + geom_segment(x= best_tag[1,]$BEG, xend =best_tag[1,]$BEG, y= best_tag[1,]$logP+0.3, yend = best_tag[1,]$logP, size = 0.01, arrow =  arrow(length = unit(0.01, "npc")), colour = "black")
  
  # save image file
  plot_out <- paste(out_dir,"ld_plot_positive_R2",trait,"_",sig$VAR,".png",sep = "")
  #print(plot_out)
  print(pp)
  ggsave(plot_out,width = 10,height = 6)
  
}


# Collect + parse command line arguments
args <- commandArgs(T)

if(length(args) < 1){
  args <- c("--help")
}

if("--help" %in% args) {
  cat("ld.snp.plot.fun.R
      Arguments:
      --snp_assoc=snp_association_file   - EPACTS format association results for the flanking SNPs
      --snp_ld=snp_ld_file   - LD between SNPs and SV genotype 
      --trait=trait   - string
      --sv_assoc=sv_association_file    - SV association results (TRAIT column required)
      --out=output_dir - the directory to put the plots
      --help              - print this text
      ")
  
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$snp_assoc) | is.null(argsL$snp_ld) | is.null(argsL$trait) | is.null(argsL$sv_assoc) | is.null(argsL$out)){
  print("missing arguements!")
  q(save = "no")
}

# setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/6-Summarize_replication/data/novel_sig_plots/40551/data/')
# assoc <- read.table('snp.ln_P_CRP_combined_rn.emmax.kin.lm.assoc.epacts.gz',header = T, comment.char = "")
assoc <- read.table(argsL$snp_assoc,header = T, comment.char = "")
# ld <- read.table('ld.snp.continous.dechr.txt',header = T)
ld <- read.table(argsL$snp_ld,header = T)
# out_dir="/Users/leichen/Desktop/Lab/Finmetseq_paper/3-Candidate_analysis/data/cand_plots/chr7_63298201_63300000/plots/"
#sig_table <- read.table('assoc.results.all_traits.sorted.txt', header = T)
sig_table <- read.table(argsL$sv_assoc, header = T)
# trait="ln_P_CRP_combined_rn"
#plot_ld_w_snps(sig_table,trait,ld,assoc,'/Users/leichen/Desktop/Lab/Finmetseq_paper/')
plot_ld_w_snps(sig_table,argsL$trait,ld,assoc,argsL$out)

