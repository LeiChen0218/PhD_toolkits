setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/highQ_anno/')
library(ggplot2)
library(tidyr)

count <- read.table('count_per_sample.cnvnator.vcf3.txt', header = T)
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

count_long <- gather(count,svtype,var_count, DUP_r2:MIXED_r2)

df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/cnvnator.sv_counts.per_sample.vcf3.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator SV counts by sample, ordered by sequencing date (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count_long <- gather(count,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

pdf('plots/cnvnator.sv_counts.per_sample.vcf3.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF3)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

