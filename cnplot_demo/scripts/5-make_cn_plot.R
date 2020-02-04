#####
# 5-make_cn_plot.R
# Given cn file, signal CNV info, signal CNV genotype, plotting copy number profile
# Format requirement for genotype file: three columns - SAMPLE,GT_DIS (sample id, discrete genotype)
# Author: Lei Chen (leichen@wustl.edu)
# Date: 2019-9-17

# load all the required packages
library(grid)
library(ggplot2)
library(gplots)
library(RColorBrewer)  
library(reshape2) 
library(gridExtra)
library(plotly)
library(ggdendro)

############
# read in tables and call functions 

setwd('/Users/leichen/Desktop/Lab/cnplot_demo/data')
cn_matrix <-  read.table('all_sample.cn',header = T)
dosage_table <- read.table('exp.20samples.list')
colnames(dosage_table) <- c("SAMPLE","GT")

get_dosage_order <- function(){
  dosage_table$SAMPLE <- as.character(dosage_table$SAMPLE)
  rep_sample_name <- gsub("^X","" ,colnames(cn_matrix))
  sample_name_clean <- gsub("[.]","-",rep_sample_name)
  rep_sample_dosage <- subset(dosage_table,dosage_table$SAMPLE %in%  sample_name_clean)
  rep_sample_dosage$SAMPLE <- factor(rep_sample_dosage$SAMPLE)
  sample_order <- rep_sample_dosage[order(rep_sample_dosage$GT),]$SAMPLE
  return(as.character(sample_order))
}

###
# Make the plot

  # select the region range to plot
  start=26766201
  end=27418201
  sub_matrix <- cn_matrix[cn_matrix$PhysicalPosition >= start & cn_matrix$PhysicalPosition <= end,]
  cn_only <- sub_matrix[,-(1:3)]
  
  # transpose and format the matrix for plotting
  matrix_t <- t(cn_only)
  colnames(matrix_t) <- sub_matrix$PhysicalPosition
  sample_name <- gsub("^X","" ,rownames(matrix_t))
  sample_name_clean <- gsub("[.]","-",sample_name)
  rownames(matrix_t) <- sample_name_clean
  
  # scale the matrix (set up the saturation color)
  cn_sat=4
  scaled_cn <- apply(matrix_t, c(1,2), FUN=function(x) min(x,cn_sat))
  
  # get the sample order for plotting
  order_ref <- get_dosage_order()
  sorted_matrix <- scaled_cn[order_ref,]
  
  # melt the matrix for plotting 
  df_sub.melted <- as.data.frame(melt(sorted_matrix))
  colnames(df_sub.melted)<-c("samples","pos","copynumber")
  df_sub.melted$samples <-as.character(df_sub.melted$samples)
  #head(df_sub.melted)
  
  # make the plot
  label_x= paste("pos on chromosome 11")
  g=ggplot(df_sub.melted,aes(x=pos,y=samples,fill=copynumber)) + geom_tile()+
    scale_fill_gradient2(low='blue', high='red',mid='white', midpoint = 2, limits= c(0,cn_sat),na.value = "red")+
    ggtitle("copy number plot demo, chr11_26969501_27010500")+xlab(label_x)+scale_y_discrete(limits = order_ref)
  
    var_pos=26969501
    var_end=27010500
    g=g+geom_vline(aes(xintercept=var_pos), colour="#990000", linetype="dashed")
    g=g+geom_vline(aes(xintercept=var_end), colour="#990000", linetype="dashed")
  print(g)

  ggsave('plots/cn_plot_demo.chr11_26969501_27010500.png',width = 14,height = 12)
  
