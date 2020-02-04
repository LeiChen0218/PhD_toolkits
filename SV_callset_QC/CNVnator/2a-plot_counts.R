setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/sample_level_qc')
library(ggplot2)
library(tidyr)

count <- read.table('count_per_sample.cnvnator.vcf1.txt', header = T)
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

count_long <- gather(count,svtype,var_count, DUP_r2:MIXED_r2)

df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/cnvnator.sv_counts.per_sample.vcf1.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count_long <- gather(count,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

pdf('plots/cnvnator.sv_counts.per_sample.vcf1.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

####
# 2- after filter out the the batch FPs

count <- read.table('count_per_sample.cnvnator.batch_fp_filtered.txt', header = T)

count_long <- gather(count,svtype,var_count, DUP_r2:MIXED_r2)

df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/cnvnator.sv_counts.per_sample.batch_fp_filtered.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count_long <- gather(count,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="nid")

pdf('plots/cnvnator.sv_counts.per_sample.batch_fp_filtered.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("CNVnator autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

###
# filter based on overall count 
count$all_r2 <- count$DEL_r2 + count$DUP_r2 + count$MIXED_r2
count$all_rm <- count$DEL_rm + count$DUP_rm + count$MIXED_rm

pdf('plots/per_sample_count.hist.pdf',width = 6, height = 4)
median_r2 = median(count$all_r2)
mad_r2 = mad(count$all_r2)
ggplot(count,aes(x=all_r2))+geom_histogram(bins=100)+ggtitle('SV count per sample, ref=2')+
  geom_vline(xintercept = c(median_r2 - mad_r2*5, median_r2 + mad_r2*5), linetype="dashed")

median_rm = median(count$all_rm)
mad_rm = mad(count$all_rm)
ggplot(count,aes(x=all_rm))+geom_histogram(bins=100)+ggtitle('SV count per sample, ref=mode')+
  geom_vline(xintercept = c(median_rm - mad_rm*5, median_rm + mad_rm*5), linetype="dashed")
dev.off()

meds <- apply(count[2:9],2,median)
mads <- apply(count[2:9],2,mad)

count_outlier <- function(v){
  med_v <- median(v)
  mad_v <- mad(v)
  outlier <- v - med_v > 5* mad_v
  return(outlier)
}

apply(count[2:9],2,count_outlier)
count$outlier <- count_outlier(count$all_rm)
summary(count$outlier)

outlier <- count[count$outlier,]

write.table(outlier, 'sample_outliers.all_rm.5mad.txt', quote = F, row.names = F,sep="\t")
