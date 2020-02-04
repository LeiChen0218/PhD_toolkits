setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/GS/sample_level_qc')
library(ggplot2)
library(tidyr)

count <- read.table('count_per_sample.gs.vcf2.txt', header = T)
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

count_long <- gather(count,svtype,var_count, DUP_r2:MIXED_r2)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

temp <- seq_date[seq_date$cram_id %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$cram_id

pdf('plots/gs.sv_counts.per_sample.vcf2.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count_long <- gather(count,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

pdf('plots/gs.sv_counts.per_sample.vcf2.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

## filtered version

count <- read.table('count_per_sample.gs.batch_fp_filtered.txt', header = T)
seq_date <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
seq_date$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',seq_date$seq.date,perl = T)

count_long <- gather(count,svtype,var_count, DUP_r2:MIXED_r2)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

temp <- seq_date[seq_date$cram_id %in% count$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$cram_id

pdf('plots/gs.sv_counts.per_sample.vcf2.batch_fp_filtered.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

count_long <- gather(count,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

pdf('plots/gs.sv_counts.per_sample.batch_fp__filtered.vcf2.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

########
# filter samples with abnormal variants count 
finns <- df[df$cohort != "Control",] 
ggplot(finns, aes(x=var_count, fill=svtype))+geom_histogram(bins=100)+
  facet_grid(.~svtype, scales = "free_x")
ggsave('plots/var_count.hist.batch_fp_filtered.png', width = 10, height = 6)

head(count)
count_finn <- count[count$sample %in% finns$sample,]
count_finn$DUP_r2_outliers <- abs(count_finn$DUP_r2 - median(count_finn$DUP_r2)) > 10*mad(count_finn$DUP_r2)
count_finn$DEL_r2_outliers <- abs(count_finn$DEL_r2 - median(count_finn$DEL_r2)) > 10*mad(count_finn$DEL_r2)
count_finn$MIXED_r2_outliers <- abs(count_finn$MIXED_r2 - median(count_finn$MIXED_r2)) > 10*mad(count_finn$MIXED_r2)
count_finn$all_r2 <- count_finn$DUP_r2+ count_finn$DEL_r2 + count_finn$MIXED_r2
count_finn$all_r2_outliers 

table(count_finn[c("DUP_r2_outliers","MIXED_r2_outliers")])

count_finn$is_outlier <- count_finn$DUP_r2_outliers | count_finn$DEL_r2_outliers  | count_finn$MIXED_r2_outliers
outliers <- merge(count_finn[count_finn$is_outlier,], seq_date, by.x = "sample", by.y = "cram_id")

table(outliers$cohort)
#Control     corogene Dyslipidemia  Eastern.Fin      FINRISK       METSIM 
#0            0           52            0           28           37

write.table(outliers, "gs.var_count.outliers.finns.table", quote = F, row.names = F,sep = "\t")
count_exclude <- count[!count$sample %in% outliers$sample,] 

count_long <- gather(count_exclude,svtype,var_count, DUP_r2:MIXED_r2)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

temp <- seq_date[seq_date$cram_id %in% count_exclude$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$cram_id

pdf('plots/gs.sv_counts.per_sample.vcf3.ref2.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()


count_long <- gather(count_exclude,svtype,var_count, DUP_rm:MIXED_rm)
df <- merge(count_long,seq_date,by.x="sample",by.y="cram_id")

temp <- seq_date[seq_date$cram_id %in% count_exclude$sample,]
sample_order <- temp[order(as.Date(temp$seq.date, format="%Y-%m-%d")),]$cram_id

pdf('plots/gs.sv_counts.per_sample.vcf3.refmode.pdf',width = 8, height = 5)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=sequencing_year))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample, ordered by sequencing date (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = sample_order)
ggplot(df,aes(x=sample,y=var_count,fill=svtype))+
  geom_bar(stat="identity")+ggtitle("GS autosome SV counts by sample (VCF2)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

