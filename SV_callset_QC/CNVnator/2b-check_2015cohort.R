setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/sample_level_qc/')
library(ggplot2)

c2015 <- read.table('cnvnator.highQ.Freq_carrier.2015.txt', header = T)
c2016 <- read.table('cnvnator.highQ.Freq_carrier.2016.txt', header = T)
c2017 <- read.table('cnvnator.highQ.Freq_carrier.2017.txt', header = T)
#head(cnvnator)

pdf('plots/subset_by_seq_year.pdf', width = 6, height = 4)
ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset))+
  geom_point(size=0.2)+ggtitle('subset: 2015')
ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset))+
  geom_point(size=0.2)+ggtitle('subset: 2016')
ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset))+
  geom_point(size=0.2)+ggtitle('subset: 2017')
dev.off()

###
# Fisher's exact test for enrichment
#[leichen@blade18-1-6 check_cohort15_artefact]$ wc -l nid.2015.list
#1316 nid.2015.list
#[leichen@blade18-1-6 check_cohort15_artefact]$ wc -l nid.2016.list
#1606 nid.2016.list
#[leichen@blade18-1-6 check_cohort15_artefact]$ wc -l nid.2017.list
#2163 nid.2017.list

size15 = 1316
size16 = 1616
size17 = 2163

c2015$ncarrier_subset <- round(c2015$carrier_rate_subset*size15)
c2015$ncarrier_nonsubset <- round(c2015$carrier_rate_nonsubset*(size16+size17))
c2015$nonc_sub <- size15 - c2015$ncarrier_subset
c2015$nonc_nonsub <-  (size16+size17 - c2015$ncarrier_nonsubset)
c2015$enrich <- c2015$carrier_rate_subset > c2015$carrier_rate_nonsubset

c2016$ncarrier_subset <- round(c2016$carrier_rate_subset*size16)
c2016$ncarrier_nonsubset <- round(c2016$carrier_rate_nonsubset*(size15+size17))
c2016$nonc_sub <- size16 - c2016$ncarrier_subset
c2016$nonc_nonsub <-  (size15+size17 - c2016$ncarrier_nonsubset)
c2016$enrich <- c2016$carrier_rate_subset > c2016$carrier_rate_nonsubset

c2017$ncarrier_subset <- round(c2017$carrier_rate_subset*size17)
c2017$ncarrier_nonsubset <- round(c2017$carrier_rate_nonsubset*(size15+size16))
c2017$nonc_sub <- size17 - c2017$ncarrier_subset
c2017$nonc_nonsub <-  (size15+size16 - c2017$ncarrier_nonsubset)
c2017$enrich <- c2017$carrier_rate_subset > c2017$carrier_rate_nonsubset

var <- c2015[1,]
run_fisher_exact <- function(var){
  table <- matrix(c(var$ncarrier_subset,var$ncarrier_nonsubset,var$nonc_sub,var$nonc_nonsub), nrow = 2)
  if(var$enrich){
    alt="greater"
  }else{
    alt="less"
  }
  test <- fisher.test(table, alternative = alt)
  #if(test$p.value==0)
  #  return(150)
  return(-log10(test$p.value))
}

c2015$logFisherP <- 1
for(i in 1:53661){
  c2015[i,]$logFisherP <- run_fisher_exact(c2015[i,])
}
#c2015$fisherP_filter <- c2015$logFisherP >= 100
ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
#ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=fisherP_filter))+
#  geom_point(size=0.2)
ggsave('plots/enrichment_test.2015.png',width = 6,height = 5)

c2015$logFisherP <- sapply(c2015$logFisherP, function(x) min(x,300))
ggplot(c2015, aes(x=logFisherP))+geom_histogram(bins=100)
ggsave('plots/enrichment_p.hist.2015.png')

c2016$logFisherP <- 1
for(i in 1:53661){
  c2016[i,]$logFisherP <- run_fisher_exact(c2016[i,])
}
#c2016$fisherP_filter <- c2016$logFisherP >= 100
ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
#ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=fisherP_filter))+
#  geom_point(size=0.2)
ggsave('plots/enrichment_test.2016.png',width = 6,height = 5)

c2016$logFisherP <- sapply(c2016$logFisherP, function(x) min(x,300))
ggplot(c2016, aes(x=logFisherP))+geom_histogram(bins=100)
ggsave('plots/enrichment_p.hist.2016.png')


c2017$logFisherP <- 1
for(i in 1:53661){
  c2017[i,]$logFisherP <- run_fisher_exact(c2017[i,])
}
#c2017$fisherP_filter <- c2017$logFisherP >= 100
ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
#ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=fisherP_filter))+
#  geom_point(size=0.2)
ggsave('plots/enrichment_test.2017.png',width = 6,height = 5)

c2017$logFisherP <- sapply(c2017$logFisherP, function(x) min(x,300))
ggplot(c2017, aes(x=logFisherP))+geom_histogram(bins=100)
ggsave('plots/enrichment_p.hist.2017.png')

all <- as.data.frame(cbind(c2015$carrier_rate_subset,
             c2015$logFisherP,
             c2016$carrier_rate_subset,
             c2016$logFisherP,
             c2017$carrier_rate_subset,
             c2017$logFisherP))
colnames(all) <- c("rate2015","neglogP2015","rate2016","neglogP2016","rate2017","neglogP2017")
all$SV <- c2015$SV_id
all$outlier_15 <- all$neglogP2015 >= 200
all$outlier_16 <- all$neglogP2016 >= 200
all$outlier_17 <- all$neglogP2017 >= 200

all$outlier <- all$outlier_15 | all$outlier_16 | all$outlier_17
table(all$outlier)
outlier <- all[all$outlier,]
c2015$outlier <- all$outlier_15
c2016$outlier <- all$outlier_16
c2017$outlier <- all$outlier_17

pdf('plots/outlier_by_NfisherP_200.pdf',width = 8, height = 6)
ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+geom_point(size=0.2)+
  ggtitle('%carriers in 2015 vs in 2016+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+geom_point(size=0.2)+
  ggtitle('%carriers in 2016 vs in 2015+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+geom_point(size=0.2)+
  ggtitle('%carriers in 2017 vs in 2016+2015 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
dev.off()

write.table(outlier[c("SV","rate2015","rate2016","rate2017")],'outlier_by_NfisherP_200.txt', quote = F, row.names = F, sep = "\t")

#########
# try t-test (not as good as fisher's exact test)
ttest15 <- read.table('cnvnator.highQ.t_test.2015.txt',header = T) 
ttest15 <- merge(ttest15,c2015,by="SV_id")
ttest15$logTP <- -log10(ttest15$p.value)

ggplot(ttest15, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logTP))+
  geom_point(size=0.2)
ggsave('plots/corr_ttest.2015.png',height = 5, width = 6)
ggplot(ttest15, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=abs(t)))+
  geom_point(size=0.2)
ggsave('plots/corr_ttest.2015.t.png',height = 5, width = 6)
