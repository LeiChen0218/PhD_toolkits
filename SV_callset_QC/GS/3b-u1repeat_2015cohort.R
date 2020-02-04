setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/GS/sample_level_qc')
library(ggplot2)
size15 = 1316
size16 = 1616
size17 = 2163

# Fisher's exact test for each cohort
# 2015
c2015 <- read.table('gs.repeat.Freq_carrier.2015.txt', header = T)

c2015$ncarrier_subset <- round(c2015$carrier_rate_subset*size15)
c2015$ncarrier_nonsubset <- round(c2015$carrier_rate_nonsubset*(size16+size17))
c2015$nonc_sub <- size15 - c2015$ncarrier_subset
c2015$nonc_nonsub <-  (size16+size17 - c2015$ncarrier_nonsubset)
c2015$enrich <- c2015$carrier_rate_subset > c2015$carrier_rate_nonsubset

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
for(i in 1:46702){
  c2015[i,]$logFisherP <- run_fisher_exact(c2015[i,])
}

ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
ggsave('plots/enrichment_test.2015.png',width = 6,height = 5)

c2015$logFisherP <- sapply(c2015$logFisherP, function(x) min(x,300))
c2015$outlier <- c2015$logFisherP >= 200
outlier <- c2015[c2015$outlier,]
ggplot(c2015, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+
  geom_point(size=0.2)+ggtitle('%carriers in 2015 vs in 2016+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggsave('plots/enrichment_test.filtered.2015.png',width = 6,height = 5)

write.table(outlier,'c2015.outlier.NfisherP_200.txt', quote = F, row.names = F, sep = "\t")

# 2016
c2016 <- read.table('gs.high_conf.Freq_carrier.2016.txt', header = T)

c2016$ncarrier_subset <- round(c2016$carrier_rate_subset*size16)
c2016$ncarrier_nonsubset <- round(c2016$carrier_rate_nonsubset*(size15+size17))
c2016$nonc_sub <- size16 - c2016$ncarrier_subset
c2016$nonc_nonsub <-  (size15+size17 - c2016$ncarrier_nonsubset)
c2016$enrich <- c2016$carrier_rate_subset > c2016$carrier_rate_nonsubset


c2016$logFisherP <- 1
for(i in 1:46702){
  c2016[i,]$logFisherP <- run_fisher_exact(c2016[i,])
}

ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
ggsave('plots/enrichment_test.2016.png',width = 6,height = 5)

c2016$logFisherP <- sapply(c2016$logFisherP, function(x) min(x,300))
c2016$outlier <- c2016$logFisherP >= 200
outlier <- c2016[c2016$outlier,]
ggplot(c2016, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+
  geom_point(size=0.2)+ggtitle('%carriers in 2016 vs in 2016+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggsave('plots/enrichment_test.filtered.2016.png',width = 6,height = 5)

write.table(outlier,'c2016.outlier.NfisherP_200.txt', quote = F, row.names = F, sep = "\t")

# 2017
c2017 <- read.table('gs.high_conf.Freq_carrier.2017.txt', header = T)

c2017$ncarrier_subset <- round(c2017$carrier_rate_subset*size17)
c2017$ncarrier_nonsubset <- round(c2017$carrier_rate_nonsubset*(size15+size16))
c2017$nonc_sub <- size17 - c2017$ncarrier_subset
c2017$nonc_nonsub <-  (size15+size16 - c2017$ncarrier_nonsubset)
c2017$enrich <- c2017$carrier_rate_subset > c2017$carrier_rate_nonsubset


c2017$logFisherP <- 1
for(i in 1:46702){
  c2017[i,]$logFisherP <- run_fisher_exact(c2017[i,])
}

ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=logFisherP))+
  geom_point(size=0.2)
ggsave('plots/enrichment_test.2017.png',width = 6,height = 5)

c2017$logFisherP <- sapply(c2017$logFisherP, function(x) min(x,300))
c2017$outlier <- c2017$logFisherP >= 200
outlier <- c2017[c2017$outlier,]
ggplot(c2017, aes(x=carrier_rate_subset, y=carrier_rate_nonsubset,color=outlier))+
  geom_point(size=0.2)+ggtitle('%carriers in 2017 vs in 2017+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggsave('plots/enrichment_test.filtered.2017.png',width = 6,height = 5)

write.table(outlier,'c2017.outlier.NfisherP_200.txt', quote = F, row.names = F, sep = "\t")


# gs <- read.table('gs.repeat.Freq_carrier.2015.txt', header = T)
# gs$in_repeat <- gs$in_u1repeat == "Yes" | gs$in_u2repeat == "Yes"
# gs$af_diff <- gs$carrier_rate_subset-gs$carrier_rate_nonsubset 
# gs$diff_sig <- gs$af_diff > 0.5
# 
# pdf('plots/gs.2015.repeat.pdf', width = 5, height = 4)
# ggplot(gs,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_u1repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('gs, % carriers, colored by intersection with u1repeat')
# ggplot(gs,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('gs, % carriers, color by intersection with u1 or u2 repeat')
# ggplot(gs,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=af_diff))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('gs, % carriers, colored by frequency differences')
# ggplot(gs,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=diff_sig))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('gs, % carriers, colored by frequency differences > 0.5')
# dev.off()
# 
# pdf('plots/diff_freq.dist.pdf')
# ggplot(gs, aes(x=af_diff))+geom_histogram(bins=100)+ggtitle('all SVs')
# dev.off()
# 
# # might just be enough to only look at u1repeats
# table(gs$diff_sig, gs$in_u1repeat)
# 
# # No   Yes
# # FALSE 38305  7731
# # TRUE      0   666
