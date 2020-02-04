setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/Lumpy/sample_level_qc')
library(ggplot2)

c2015 <- read.table('lumpy.mie.pass.autosome.Freq_carrier.2015.txt', header = T)
size15 = 1316
size16 = 1616
size17 = 2163

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
for(i in 1:43000){
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
c2016 <- read.table('lumpy.mie.pass.autosome.Freq_carrier.2016.txt', header = T)

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
  geom_point(size=0.2)+ggtitle('%carriers in 2016 vs in 2015+2017 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggsave('plots/enrichment_test.filtered.2016.png',width = 6,height = 5)

write.table(outlier,'c2016.outlier.NfisherP_200.txt', quote = F, row.names = F, sep = "\t")

# 2017
c2017 <- read.table('lumpy.mie.pass.autosome.Freq_carrier.2017.txt', header = T)

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
  geom_point(size=0.2)+ggtitle('%carriers in 2017 vs in 2016+2015 (seq.date)\nFilter: -logP(Fisher.Exact.P) >= 200')
ggsave('plots/enrichment_test.filtered.2017.png',width = 6,height = 5)

write.table(outlier,'c2017.outlier.NfisherP_200.txt', quote = F, row.names = F, sep = "\t")



# bnd <- read.table('rm_repeat/bnd.repeat.Freq_carrier.2015.txt', header = T)
# bnd$in_repeat <- bnd$in_u1repeat == "Yes" | bnd$in_u2repeat == "Yes"
# bnd$af_diff <- bnd$carrier_rate_subset-bnd$carrier_rate_nonsubset 
# bnd$diff_sig <- bnd$af_diff > 0.5
# 
# pdf('plots/bnd.2015.repeat.pdf', width = 5, height = 4)
# ggplot(bnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_u1repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('BND, % carriers, colored by intersection with u1repeat')
# ggplot(bnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('BND, % carriers, color by intersection with u1 or u2 repeat')
# ggplot(bnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=af_diff))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('BND, % carriers, colored by frequency differences')
# ggplot(bnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=diff_sig))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('BND, % carriers, colored by frequency differences > 0.5')
# dev.off()
# 
# nonbnd <- read.table('rm_repeat/nonbnd.repeat.Freq_carrier.2015.txt', header = T)
# nonbnd$in_repeat <- nonbnd$in_u1repeat == "Yes" | nonbnd$in_u2repeat == "Yes"
# nonbnd$af_diff <- nonbnd$carrier_rate_subset-nonbnd$carrier_rate_nonsubset 
# nonbnd$diff_sig <- nonbnd$af_diff > 0.5
# 
# pdf('plots/nonbnd.2015.repeat.pdf', width = 5, height = 4)
# ggplot(nonbnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_u1repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('nonBND % carriers, colored by intersection with u1repeat')
# ggplot(nonbnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=in_repeat))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('nonBND % carriers, color by intersection with u1 or u2 repeat')
# ggplot(nonbnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=af_diff))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('nonBND % carriers, colored by frequency differences')
# ggplot(nonbnd,aes(x=carrier_rate_subset, y=carrier_rate_nonsubset, color=diff_sig))+
#   geom_point(size=0.3)+xlab('2015 cohort')+ylab('2016+2017 cohort')+
#   ggtitle('nonBND % carriers, colored by frequency differences > 0.5')
# dev.off()
# 
# all_freq <- read.table('rm_repeat/lumpy.mie.pass.autosome.Freq_carrier.2015.txt', header = T)
# all_freq$freq_diff <- all_freq$carrier_rate_subset - all_freq$carrier_rate_nonsubset
# pdf('plots/diff_freq.dist.pdf')
# ggplot(all_freq, aes(x=freq_diff))+geom_histogram(bins=100)+ggtitle('all SVs')
# ggplot(bnd, aes(x=af_diff))+geom_histogram(bins=100)+ggtitle('bnd')
# ggplot(nonbnd, aes(x=af_diff))+geom_histogram(bins=100)+ggtitle('nonbnd')
# dev.off()
# 
# # might just be enough to only look at u1repeats
# table(bnd$diff_sig, bnd$in_u1repeat)
# 
# #No  Yes
# #FALSE 9807  380
# #TRUE   233  345
# table(bnd$diff_sig, bnd$in_repeat)
# 
# #FALSE TRUE
# #FALSE  9197  990
# #TRUE    223  355
