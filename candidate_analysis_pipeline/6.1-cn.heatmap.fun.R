library(grid)
library(ggplot2)
library(gplots)
library(RColorBrewer)  
library(reshape2) 
library(gridExtra)
library(plotly)
library(ggdendro)

# setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/3-Candidate_analysis/data/cand_plots/chr2_1343601_1344800/data')
# cn_matrix <- read.table('rep100sample.cn',header = T)
# sig <- read.table('assoc.results.all_traits.sorted.txt', header = T)
# sig_pos <- 1343601        
# sig_end <- 1344800
# start <- 1337651
# end <- 1350750
# dosage_table <- read.table('genotype.t.list',header = T)  

get_clustered_order <- function(df){
  cn_dist <- dist(df, method = "euclidean")
  h.fit <- hclust(cn_dist, method="ward")
  reorder <- h.fit$order
  return(rownames(df[reorder,])) 
}

get_dosage_order <- function(cn_matrix,dosage_table){
  dosage_table$SAMPLE <- as.character(dosage_table$SAMPLE)
  rep_sample_name <- gsub("^X","" ,colnames(cn_matrix))
  sample_name_clean <- gsub("[.]","-",rep_sample_name)
  rep_sample_dosage <- subset(dosage_table,dosage_table$SAMPLE %in%  sample_name_clean)
  #tmp <- subset(sample_name_clean, !(sample_name_clean %in% dosage_table$SAMPLE))
  rep_sample_dosage$SAMPLE <- factor(rep_sample_dosage$SAMPLE)
  sample_order <- rep_sample_dosage[order(rep_sample_dosage$GT_DIS),]$SAMPLE
  return(as.character(sample_order))
}

cn_heatmap_fun <- function(cn_matrix,
                           order_ref=NULL,
                           start,
                           end,
                           cn_sat=NULL,
                           sig,
                           out){
  if(is.null(cn_sat))
    cn_sat = 4
  else
    cn_sat = as.numeric(cn_sat)
  
  sub_matrix <- cn_matrix[cn_matrix$PhysicalPosition >= start & cn_matrix$PhysicalPosition <= end,]
  
  cn_only <- sub_matrix[,-(1:3)]
  
  # transpose and format the matrix for plotting
  matrix_t <- t(cn_only)
  colnames(matrix_t) <- sub_matrix$PhysicalPosition
  sample_name <- gsub("^X","" ,rownames(matrix_t))
  sample_name_clean <- gsub("[.]","-",sample_name)
  rownames(matrix_t) <- sample_name_clean
  
  # scale the matrix (set up the saturation color) 
  scaled_cn <- apply(matrix_t, c(1,2), FUN=function(x) min(x,cn_sat))
  
  # get the sample order for plotting
  title=paste("copy number heatmap (sample ordered by SV genotype)\nsample size=100, cn saturated at",cn_sat,")",sep = " ")
  if(is.null(order_ref)){
    order_ref <- get_clustered_order(scaled_cn)
    title=paste("copy number heatmap (sample ordered by clustering)\nsample size=100, cn saturated at",cn_sat,")",sep = " ")
  }
  sorted_matrix <- scaled_cn[order_ref,]
  
  # melt the matrix for plotting 
  df_sub.melted <- as.data.frame(melt(sorted_matrix))
  colnames(df_sub.melted)<-c("samples","pos","copynumber")
  df_sub.melted$samples <-as.character(df_sub.melted$samples)
  #head(df_sub.melted)
  
  var.info=paste(sig$VAR,sig$TYPE,"length(bp)=",sig$LEN,"AF=",sig$AF)
  title=paste(var.info,title,sep="\n")
  # make the plot
  label_x= paste("pos on chromosome",sig$CHR)
  g=ggplot(df_sub.melted,aes(x=pos,y=samples,fill=copynumber)) + geom_tile()+
    scale_fill_gradient2(low='blue', high='red',mid='white', midpoint = 2, limits= c(0,cn_sat),na.value = "red")+
    ggtitle(title)+xlab(label_x)+scale_y_discrete(limits = order_ref)
  
  if(!is.null(sig))
    g=g+geom_vline(aes(xintercept=sig$POS), colour="#990000", linetype="dashed")
  if(!is.null(sig))
    g=g+geom_vline(aes(xintercept=sig$END), colour="#990000", linetype="dashed")
  print(g)
  #out="cn.100bp.44368.rep100samp.27705000_27710000"
  ggsave(paste(out,".png",sep = ""),width = 14,height = 12)
  
}