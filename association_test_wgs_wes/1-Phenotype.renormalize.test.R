setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/2-Association_test/data/general/')
library(ggplot2)

pheno <- read.table('../../../general_info/qt_finnseq_20161017.5k_ids.ped', header = T, comment.char = "")
s_table <- read.table('../../../general_info/finn.samp.all.wcramdir.seqdate.table', header = T)
traits <- read.table('../../../general_info/wes.qt64.list')
sample_count <- read.table('../../../1-Callset_QC/data/general/spacial_cluster/ClusterSampleCount.lumpy.cnv.txt', header = T)
sample_count <- merge(sample_count,s_table,by.x = "SAMPLE", by.y = "cram_id")
sub_pheno <- pheno[pheno$IND_ID %in% sample_count$nid,]

#all_wes <- read.table('../../../general_info/qt_finnseq_20161017.ped', header = T, comment.char = "")
#all_wes <- merge(all_wes, s_table, by.x = "IND_ID", by.y = "oid")
#sum(all_wes$nid %in% sample_count$SAMPLE)

rank.based.INT <- function(v){
  r = rank(v)
  c=3/8 #parameter to prevent max/min from being converted to +/- infinite value 
  rank_score <- (r-c)/(length(r)-2*c+1)
  return(qnorm(rank_score)) 
}

#
re_normalize_trait <- function(trait){
  sub <- sub_pheno[c("IND_ID",trait)]
  sub <- sub[complete.cases(sub),]
  sub$renomalized <- rank.based.INT(sub[trait])
  cname <- paste(trait,"rn",sep = "_")
  colnames(sub) <- c("IND_ID",trait,cname)
  return(sub)
}


header_cols <- sub_pheno[1:22] 
trait_list <-traits$V1

for(trait in trait_list){
  sub <- re_normalize_trait(trait)
  header_cols <- merge(header_cols, sub, by = "IND_ID", all.x = T)
  #print(dim(header_cols))
}

sample_size = apply(header_cols, 2, function(x) sum(!is.na(x)))

write.table(header_cols, "qt_finnseq_20190620.renormalized.wgs.5k_ids.ped",
            quote = F, row.names = F, sep="\t")

###
# phenotype redundancy
sub <- header_cols[as.vector(trait_list)]
corr_mat <- cor(sub, use="pairwise.complete.obs")

write.table(corr_mat,'pheno.corr.mat.txt', quote = F, sep = "\t", row.names = F, col.names = F)

## renormalize all traits passed power > 0.8 (N>374.6)
samp_size <- apply(sub_pheno,2,function(x) sum(!is.na(x)))
power_pheno <- sub_pheno[,samp_size > 374]
  
header_cols <- power_pheno[1:22] 
trait_list <-colnames(power_pheno)[23:138]


for(trait in trait_list){
  sub <- re_normalize_trait(trait)
  header_cols <- merge(header_cols, sub, by = "IND_ID", all.x = T)
  #print(dim(header_cols))
}

sample_size = apply(header_cols, 2, function(x) sum(!is.na(x)))

write.table(header_cols, "qt_finnseq_20190620.renormalized.wgs.5k_ids.power80.ped",
            quote = F, row.names = F, sep="\t")

###
# phenotype redundancy
rn_traits <- grep("_rn",colnames(header_cols))
sub <- header_cols[rn_traits]
corr_mat <- cor(sub, use="pairwise.complete.obs")
apply(corr_mat, 1, function(x) sum(is.na(x)))
corr_mat[is.na(corr_mat)] <- 1

write.table(corr_mat,'pheno.corr.mat.power80.txt', quote = F, sep = "\t", row.names = F, col.names = F)


