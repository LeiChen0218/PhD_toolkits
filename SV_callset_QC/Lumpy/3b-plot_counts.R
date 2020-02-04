setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/Lumpy/sample_level_qc')
library(ggplot2)
library(tidyr)

count <- read.table('var_count.per_sample.vcf2.txt')
colnames(count) <- c("sample","svtype","homref","var_count")
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

df <- merge(count,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/lumpy.sv_counts.per_sample.vcf2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count <- read.table('var_count.per_sample.vcf2.autsome.txt')
colnames(count) <- c("sample","svtype","homref","var_count")
df <- merge(count,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/lumpy.sv_counts.per_sample.vcf2.autosome.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

# after site filtering 

count <- read.table('var_count.per_sample.vcf2.batch_fp_filtered.txt')
colnames(count) <- c("sample","svtype","homref","var_count")

df <- merge(count,seq_date,by.x="sample",by.y="nid")

temp <- seq_date[seq_date$nid %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$nid

pdf('plots/lumpy.sv_counts.per_sample.batch_fp_filtered.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("LUMPY SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

##############
# filter sample with outlier 
finns <- df[df$cohort != "Control",]
ggplot(finns, aes(x=var_count, fill=svtype))+geom_histogram(bins=100)+
  facet_grid(.~svtype, scales = "free_x")
ggsave('plots/finns_only.var_count.by_type.hist.batch_fp_filtered.png', width = 12, height = 6)

types=unique(df$svtype)
for(t in types){
  sub <- finns[finns$svtype==t,]
  sub$diff <- abs(sub$var_count - median(sub$var_count))
  sub_mad <- mad(sub$var_count)
  outlier <- sub[sub$diff > sub_mad*10,]
  print(outlier)
}
t="BND"
sub <- finns[finns$svtype==t,]
sub$diff <- abs(sub$var_count - median(sub$var_count))
sub_mad <- mad(sub$var_count)
outlier <- sub[sub$diff > sub_mad*10,]
# write three lumpy outlier samples
write.table(outlier,'lumpy.var_count.outliers.finns.table', quote=F, row.names=F,sep='\t')
