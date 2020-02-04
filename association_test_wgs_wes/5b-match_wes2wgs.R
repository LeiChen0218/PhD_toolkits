setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/3-Candidate_analysis/data/wes_match/')
library(ggplot2)

# read in wes results
bam1 <- read.table('cands.traits.bam1.txt', header = T)
bam2 <- read.table('cands.traits.bam2.txt', header = T)
meta <- read.table('cands.traits.r01_exon.meta.txt', header = T)

# read wes to cnv chain file
chain1 <- read.table('targets/cnv.exon.r_01.bams1.table', header = F)
chain2 <- read.table('targets/cnv.exon.r_01.bams2.table', header = F)
colnames(chain1) <- c("ID","REGION","R2")
colnames(chain2) <- c("ID","REGION","R2")

# read wgs results
gs <- read.table('../cand_p3/gs.candidate.p3.txt', header = T)
lumpy <- read.table('../cand_p3/lumpy.candidate.p3.txt', header = T)
cnvnator <- read.table('../cand_p3/cnvnator.candidate.p3.txt', header = T)
wgs <- as.data.frame(rbind(gs, lumpy, cnvnator))

# merge wes results
#wes <- merge(bam1,bam2, by = c("TRAIT","REGION"), all = T)
#wes <- merge(wes, meta, by = c("TRAIT","REGION"), all = T)
#wes_sub <- wes[c(1,2,4,5,7,8,10,11,13,14)]

# match wgs candidate and exons
wgs1 <- merge(chain1, wgs, by = "ID")
colnames(wgs1) <- c("CNV","REGION","R2","TRAIT_RN","CHR","POS","WGS_P","WGS_BETA","AC","AF","N")
wgs1$TRAIT <- gsub('_rn','',wgs1$TRAIT_RN)

# merge wes and wgs results
combined_bam1 <- merge(bam1,wgs1, by =  c("TRAIT","REGION"))
#write.table(all,'wes_wgs_matched.results.txt', row.names = F, sep = '\t', quote = F)

library('metap')
combined_bam1$ED <- (combined_bam1$BETA*combined_bam1$WGS_BETA) > 0 
combined_bam1$CP <- 1

for(i in 1:dim(combined_bam1)[1]){
  combined_bam1[i,]$CP <- sumlog(c(combined_bam1[i,]$PVALUE, combined_bam1[i,]$WGS_P))$p
}

# bams2
wgs2 <- merge(chain2, wgs, by = "ID")
colnames(wgs2) <- c("CNV","REGION","R2","TRAIT_RN","CHR","POS","WGS_P","WGS_BETA","AC","AF","N")
wgs2$TRAIT <- gsub('_rn','',wgs2$TRAIT_RN)
combined_bam2 <- merge(bam2,wgs2, by =  c("TRAIT","REGION"))
combined_bam2$ED <- (combined_bam2$BETA*combined_bam2$WGS_BETA) > 0 
combined_bam2$CP <- 1
for(i in 1:dim(combined_bam2)[1]){
  combined_bam2[i,]$CP <- sumlog(c(combined_bam2[i,]$PVALUE, combined_bam2[i,]$WGS_P))$p
}

# meta
chain <-read.table('targets/cnv.exon.r_01.table', header = F)
colnames(chain) <- c("ID","REGION")
wgs_meta <- merge(chain, wgs, by = "ID")
colnames(wgs_meta) <- c("CNV","REGION","TRAIT_RN","CHR","POS","WGS_P","WGS_BETA","AC","AF","N")
wgs_meta$TRAIT <- gsub('_rn','',wgs_meta$TRAIT_RN)
combined_meta <- merge(meta,wgs_meta, by =  c("TRAIT","REGION"))
combined_meta$ED_re <- (combined_meta$BETA_RE*combined_meta$WGS_BETA) > 0 
combined_meta$CP_re <- 1
combined_meta$ED_fe <- (combined_meta$BETA_FE*combined_meta$WGS_BETA) > 0 
combined_meta$CP_fe <- 1
for(i in 1:dim(combined_meta)[1]){
  combined_meta[i,]$CP_re <- sumlog(c(combined_meta[i,]$PVALUE_RE, combined_meta[i,]$WGS_P))$p
  combined_meta[i,]$CP_fe <- sumlog(c(combined_meta[i,]$PVALUE_FE, combined_meta[i,]$WGS_P))$p
}


ggplot(combined_bam1, aes(x=-log10(CP)))+geom_histogram(bins=100)+
  ggtitle('combined p distribution, wes batch1 ')
ggplot(combined_bam2, aes(x=-log10(CP)))+geom_histogram(bins=100)+
  ggtitle('combined p distribution, wes batch2 ')
ggplot(combined_meta, aes(x=-log10(CP_fe)))+geom_histogram(bins=100)+
  ggtitle('combined p distribution, meta, fixed effect ')
ggplot(combined_meta, aes(x=-log10(CP_re)))+geom_histogram(bins=100)+
  ggtitle('combined p distribution, meta, random effect ')

write.table(combined_bam1, 'wes_wgs_matched.bam1.txt',row.names = F, sep = '\t', quote = F)
write.table(combined_bam2, 'wes_wgs_matched.bam2.txt',row.names = F, sep = '\t', quote = F)
write.table(combined_meta, 'wes_wgs_matched.meta.txt',row.names = F, sep = '\t', quote = F)

valid1 <- combined_bam1[combined_bam1$ED & combined_bam1$CP < 0.00000189,]
valid2 <- combined_bam2[combined_bam2$ED & combined_bam2$CP < 0.00000189,]
combined_meta$valid_fe <- combined_meta$ED_fe & combined_meta$CP_fe < 0.00000189 
combined_meta$valid_re <- combined_meta$ED_re & combined_meta$CP_re < 0.00000189 
valid_meta <- combined_meta[combined_meta$valid_fe,] 

write.table(valid1,'wes_wgs_matched.valid.bam1.txt', row.names = F, sep = '\t', quote = F)
write.table(valid2,'wes_wgs_matched.valid.bam2.txt', row.names = F, sep = '\t', quote = F)
write.table(valid_meta,'wes_wgs_matched.valid.meta.txt', row.names = F, sep = '\t', quote = F)

valid1 <- combined_bam1[combined_bam1$ED & combined_bam1$CP < 0.00001,]
valid2 <- combined_bam2[combined_bam2$ED & combined_bam2$CP < 0.00001,]
combined_meta$valid_fe <- combined_meta$ED_fe & combined_meta$CP_fe < 0.00001 
combined_meta$valid_re <- combined_meta$ED_re & combined_meta$CP_re < 0.00001
valid_meta <- combined_meta[combined_meta$valid_fe,] 

write.table(valid1,'wes_wgs_matched.subthre.bam1.txt', row.names = F, sep = '\t', quote = F)
write.table(valid2,'wes_wgs_matched.subthre.bam2.txt', row.names = F, sep = '\t', quote = F)
write.table(valid_meta,'wes_wgs_matched.subthre.meta.txt', row.names = F, sep = '\t', quote = F)
