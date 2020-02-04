setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/IRS_tuning/')
library(ggplot2)
library(tidyr)
source("/Users/leichen/Documents/Scripts box/R/multiplot.R")

data <- read.table('rare.tuning.table', header = T)

get_fdr <- function(df, min_prob = 1){
  sub = df[!is.na(df$IRS_P),]
  sub = sub[sub$NPROBES >= min_prob,]
  denominator = dim(sub)[1]
  nominator <- dim(sub[sub$IRS_P >= 0.5,])[1]
  return(min(c(nominator*2/denominator),1))
}

#temp <- data[!is.na(data$IRS_P),]
#temp$fp <- temp$IRS_P >= 0.5
#get_fdr(temp[temp$NNONREF > 10,], min_prob = 2)
#ggplot(temp[temp$NNONREF > 10,], aes(x=DIST))+geom_histogram(bins=100)
#ggplot(temp, aes(x=fp,y=NNONREF))+geom_violin()+ylim(0,20)+facet_grid(.~SVTYPE)

filter_callset <- function(var_type="All",
                           freq_thred= 0) {
  if(var_type=="All")
    sub <- data[  data$NNONREF >= freq_thred,]
  else
    sub <- data[  data$SVTYPE == var_type & data$NNONREF >= freq_thred,]
  return(sub)
}

get_overlap_frac <- function(df, callset = "IN_SV1KG"){
  sub = df[callset]
  rm_na = sub[!is.na(sub)]
  denominator <- length(rm_na)
  nominator <- sum(rm_na == "Yes")
  return(nominator/denominator)
}

tune_by_ncarrier_dup <- data.frame(
  freq_thred =seq(1,100,by=1),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

tune_by_ncarrier_del <- data.frame(
  freq_thred =seq(1,100,by=1),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

tune_by_ncarrier_mixed <- data.frame(
  freq_thred =seq(1,100,by=1),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

for(i in 1:100){
  filtered_df <- filter_callset(freq_thred = i, var_type = "DUP")
  tune_by_ncarrier_dup[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_ncarrier_dup[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_ncarrier_dup[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
  filtered_df <- filter_callset(freq_thred = i, var_type = "DEL")
  tune_by_ncarrier_del[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_ncarrier_del[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_ncarrier_del[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
  filtered_df <- filter_callset(freq_thred = i, var_type = "MIXED")
  tune_by_ncarrier_mixed[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_ncarrier_mixed[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_ncarrier_mixed[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
}

tune_by_ncarrier_dup_long <- gather(tune_by_ncarrier_dup, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)
tune_by_ncarrier_mixed_long <- gather(tune_by_ncarrier_mixed, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)
tune_by_ncarrier_del_long <- gather(tune_by_ncarrier_del, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)

pdf('plots/rare.tune_by_freq.by_sv_type.pdf', width = 8, height = 6)
ggplot(tune_by_ncarrier_dup_long,aes(x=freq_thred, y=fraction, color=measure))+geom_point()+ggtitle("DUP")+xlim(0,20)
ggplot(tune_by_ncarrier_del_long,aes(x=freq_thred, y=fraction, color=measure))+geom_point()+ggtitle("DEL")+xlim(0,20)
ggplot(tune_by_ncarrier_mixed_long,aes(x=freq_thred, y=fraction, color=measure))+geom_point()+ggtitle("MIXED")+xlim(0,20)
dev.off()

get_fdr(filter_callset(freq_thred = 2), min_prob = 2)

dup<-filter_callset(freq_thred = 2,var_type = "DUP") # FDR < 0.1
del<-filter_callset(freq_thred = 5,var_type = "DEL") # min FDR
mixed<-filter_callset(freq_thred = 7,var_type = "MIXED") # FDR < 0.1

all <- as.data.frame(rbind(dup,del,mixed))
get_fdr(all, min_prob = 2)
# 0.09355976
