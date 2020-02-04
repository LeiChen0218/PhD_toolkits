setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/general/spacial_cluster/')

gs_count <- read.table('ClusterSampleCount.gs.cnv.txt', header = T)
cnvnator_count <- read.table('ClusterSampleCount.cnvnator.txt', header = T)
lumpy_count <- read.table('ClusterSampleCount.lumpy.cnv.txt', header = T)
samp_info <- read.table('../../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
samp_info$sequencing_year <-  gsub('-[0-9]+-[0-9]+','',samp_info$seq.date,perl = T)

library(tidyr)
gs_df <- merge(gs_count, samp_info, by.x = "SAMPLE", by.y = "nid")
cnvnator_df <- merge(cnvnator_count, samp_info, by.x = "SAMPLE", by.y = "nid")
lumpy_df <- merge(lumpy_count, samp_info,  by.x = "SAMPLE", by.y = "cram_id")
library(ggplot2)

cnvnator_order <- cnvnator_df[order(as.Date(cnvnator_df$seq.date, format="%Y-%m-%d")),]$SAMPLE

pdf('plots/cnvnator.cluster.per_sample.highQ.pdf',width = 8, height = 5)
ggplot(cnvnator_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (CNVNATOR), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = cnvnator_order)
ggplot(cnvnator_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (CNVNATOR)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
ggplot(cnvnator_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (CNVNATOR), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = cnvnator_order)
ggplot(cnvnator_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (4936 Finns)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()


gs_order <- gs_df[order(as.Date(gs_df$seq.date, format="%Y-%m-%d")),]$SAMPLE

pdf('plots/gs.cluster.per_sample.highQ.pdf',width = 8, height = 5)
ggplot(gs_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (GS), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = gs_order)
ggplot(gs_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (GS)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
ggplot(gs_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (GS), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = gs_order)
ggplot(gs_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (GS)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()

lumpy_order <- lumpy_df[order(as.Date(lumpy_df$seq.date, format="%Y-%m-%d")),]$SAMPLE

pdf('plots/lumpy.cluster.per_sample.highQ.pdf',width = 8, height = 5)
ggplot(lumpy_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (lumpy), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = lumpy_order)
ggplot(lumpy_df,aes(x=SAMPLE,y=COUNT_NONREF,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONREF Cluster counts by sample (lumpy)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
ggplot(lumpy_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=sequencing_year ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (lumpy), ordered by sequencing date")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_x_discrete(limits = lumpy_order)
ggplot(lumpy_df,aes(x=SAMPLE,y=COUNT_NONVAR,fill=cohort ))+
  geom_bar(stat="identity")+ggtitle("NONVAR Cluster counts by sample (lumpy)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_grid(.~cohort, scales = "free_x")
dev.off()



