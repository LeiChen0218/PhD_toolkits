setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/general/spacial_cluster')

library(ggplot2)

gs <- read.table('gs.spacial_cluster.autosome.var_count.var_type_count.txt')
cnvnator <- read.table('cnvnator.spacial_cluster.var_count.var_type_count.txt')
lumpy <- read.table('lumpy.spacial_cluster.var_count.var_type_count.txt')
lumpy_cnv <- read.table('lumpy.cnv.spacial_cluster.var_count.var_type_count.txt')
colnames(gs) <- c("cluster","var_count","var_type_count")
colnames(cnvnator) <- c("cluster","var_count","var_type_count")
colnames(lumpy) <- c("cluster","var_count","var_type_count")
colnames(lumpy_cnv) <- c("cluster","var_count","var_type_count")
head(gs)

gs$callset <- "gs"
lumpy$callset <- "lumpy"
lumpy_cnv$callset <- "lumpy(del+dup)"
cnvnator$callset <- "cnvnator"

mean(gs$var_count)
mean(cnvnator$var_count)
mean(lumpy$var_count)
mean(lumpy_cnv$var_count)

max(gs$var_count)
max(cnvnator$var_count)
max(lumpy$var_count)
max(lumpy_cnv$var_count)
 
sum(lumpy$var_count == 1)/length(lumpy$cluster)
sum(lumpy_cnv$var_count == 1)/length(lumpy_cnv$cluster)
sum(gs$var_count == 1)/length(gs$cluster)
sum(cnvnator$var_count == 1)/length(cnvnator$cluster)

df <- as.data.frame(rbind(gs,lumpy,lumpy_cnv,cnvnator))
df$var_count_merged <- sapply(df$var_count, function(x) min(x,10))

ggplot(df, aes(x=var_count_merged, fill=callset))+geom_histogram(bins=10)+
  facet_grid(callset~.)+ggtitle('distribution of #. var in each cluster, saturate at 10')
ggsave('plots/n_var.dist.per.d10cluster.png',width = 6, height = 4)

