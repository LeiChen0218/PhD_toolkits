setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/Lumpy/highQ_anno/')
library(ggplot2)
library(tidyr)

count <- read.table('var_count.per_sample.autosome.vcf4.txt')
colnames(count) <- c("sample","svtype","homref","var_count")
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

df <- merge(count,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/lumpy.sv_counts.per_sample.autosome.vcf4.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample, ordered by sequencing date (VCF4)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample (VCF4)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

