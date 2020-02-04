setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/GS')
library(ggplot2)
library(tidyr)
source("/Users/leichen/Documents/Scripts box/R/multiplot.R")

gs <- read.table("tuning.table", header = T)

filter_callset <- function(df,
                           var_type="All",
                           qual_score_thred = 0) {
  
  if(var_type=="All")
    sub <- df[  df$qual >= qual_score_thred,]
  else
      sub <- df[  df$SVTYPE == var_type & df$QUAL >= qual_score_thred,]
  
  return(sub)
}

get_fdr <- function(df, min_prob = 1){
  sub = df[!is.na(df$IRS_P),]
  sub = sub[sub$NPROBES >= min_prob,]
  denominator = dim(sub)[1]
  nominator <- dim(sub[sub$IRS_P >= 0.5,])[1]
  return(min(c(nominator*2/denominator),1))
}

get_overlap_frac <- function(df, callset = "in_gs5k"){
  sub = df[callset]
  rm_na = sub[!is.na(sub)]
  denominator <- length(rm_na)
  nominator <- sum(rm_na == "Yes")
  return(nominator/denominator)
}


#######
# gs 
# by sv type
tune_by_gscnqual_dup <- data.frame(
  qual_score_thred = seq(0,9.9,by=0.1),
  fdr_all = rep(1,100),
  in_sv1kg = rep(0,100),
  in_gnomad = rep(0,100)
)

tune_by_gscnqual_del <- data.frame(
  qual_score_thred = seq(0,9.9,by=0.1),
  fdr_all = rep(1,100),
  in_sv1kg = rep(0,100),
  in_gnomad = rep(0,100)
)

tune_by_gscnqual_mixed <- data.frame(
  qual_score_thred = seq(0,9.9,by=0.1),
  fdr_all = rep(1,100),
  in_sv1kg = rep(0,100),
  in_gnomad = rep(0,100)
)

for(i in 1:100){
  filtered_df <- filter_callset(gs,qual_score_thred = (i-1)*0.1, var_type = "DUP")
  tune_by_gscnqual_dup[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_gscnqual_dup[i,]$in_sv1kg <- get_overlap_frac(filtered_df,"IN_SV1KG" )
  tune_by_gscnqual_dup[i,]$in_gnomad <- get_overlap_frac(filtered_df,"IN_GNOMAD" )
  filtered_df <- filter_callset(gs,qual_score_thred = (i-1)*0.1, var_type = "DEL")
  tune_by_gscnqual_del[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_gscnqual_del[i,]$in_sv1kg <- get_overlap_frac(filtered_df,"IN_SV1KG" )
  tune_by_gscnqual_del[i,]$in_gnomad <- get_overlap_frac(filtered_df,"IN_GNOMAD" )
  filtered_df <- filter_callset(gs,qual_score_thred = (i-1)*0.1, var_type = "MIXED")
  tune_by_gscnqual_mixed[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_gscnqual_mixed[i,]$in_sv1kg <- get_overlap_frac(filtered_df,"IN_SV1KG" )
  tune_by_gscnqual_mixed[i,]$in_gnomad <- get_overlap_frac(filtered_df,"IN_GNOMAD" )
}

tune_by_gscnqual_dup_long <- gather(tune_by_gscnqual_dup, measure, fraction, fdr_all:in_gnomad, factor_key=TRUE)
tune_by_gscnqual_mixed_long <- gather(tune_by_gscnqual_mixed, measure, fraction, fdr_all:in_gnomad, factor_key=TRUE)
tune_by_gscnqual_del_long <- gather(tune_by_gscnqual_del, measure, fraction, fdr_all:in_gnomad, factor_key=TRUE)

pdf('plots/tune_by_qual_score.by_sv_type.pdf', width = 8, height = 6)
ggplot(tune_by_gscnqual_dup_long,aes(x=qual_score_thred, y=fraction, color=measure))+geom_point()+ggtitle("DUP")
ggplot(tune_by_gscnqual_del_long,aes(x=qual_score_thred, y=fraction, color=measure))+geom_point()+ggtitle("DEL")
ggplot(tune_by_gscnqual_mixed_long,aes(x=qual_score_thred, y=fraction, color=measure))+geom_point()+ggtitle("MIXED")
dev.off()




