setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/IRS_tuning/')
library(ggplot2)
library(tidyr)
source("/Users/leichen/Documents/Scripts box/R/multiplot.R")

data <- read.table('common.tuning.table', header = T)

filter_callset <- function(var_type="All",
                           MEAN_SEP_thred= 0,
                           DIP_P_thred = 1) {
  if(var_type=="All")
    sub <- data[  data$DIP_P <= DIP_P_thred & data$MEAN_SEP  >= MEAN_SEP_thred,]
  else
    sub <- data[  data$SVTYPE == var_type & data$DIP_P <= DIP_P_thred & data$MEAN_SEP  >= MEAN_SEP_thred,]
  return(sub)
}

get_fdr <- function(df, min_prob = 1){
  sub = df[!is.na(df$IRS_P),]
  sub = sub[sub$NPROBES >= min_prob,]
  denominator = dim(sub)[1]
  nominator <- dim(sub[sub$IRS_P >= 0.5,])[1]
  return(min(c(nominator*2/denominator),1))
}

get_overlap_frac <- function(df, callset = "IN_SV1KG"){
  sub = df[callset]
  rm_na = sub[!is.na(sub)]
  denominator <- length(rm_na)
  nominator <- sum(rm_na == "Yes")
  return(nominator/denominator)
}

# dip p? -- probably no need to tune on it separately 
# temp <- data[data$DIP_P>0.5 & data$MEAN_SEP>0.47,]
# get_fdr(temp)
#summary(temp$IRS_P)

tune_by_meansep_dup <- data.frame(
  mean_sep_thred =seq(0.01,1,by=0.01),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

tune_by_meansep_del <- data.frame(
  mean_sep_thred =seq(0.01,1,by=0.01),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

tune_by_meansep_mixed <- data.frame(
  mean_sep_thred =seq(0.01,1,by=0.01),
  fdr_all = rep(1,100),
  frac_SV1kg = rep(0,100),
  frac_gnomad = rep(0,100)
)

for(i in 1:100){
  filtered_df <- filter_callset(MEAN_SEP_thred = i*0.01, var_type = "DUP")
  tune_by_meansep_dup[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_meansep_dup[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_meansep_dup[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
  filtered_df <- filter_callset(MEAN_SEP_thred = i*0.01, var_type = "DEL")
  tune_by_meansep_del[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_meansep_del[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_meansep_del[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
  filtered_df <- filter_callset(MEAN_SEP_thred = i*0.01, var_type = "MIXED")
  tune_by_meansep_mixed[i,]$fdr_all <- get_fdr(filtered_df, min_prob=2)
  tune_by_meansep_mixed[i,]$frac_SV1kg <- get_overlap_frac(filtered_df,callset = "IN_SV1KG")
  tune_by_meansep_mixed[i,]$frac_gnomad <- get_overlap_frac(filtered_df,callset = "IN_GNOMAD")
}

tune_by_meansep_dup_long <- gather(tune_by_meansep_dup, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)
tune_by_meansep_mixed_long <- gather(tune_by_meansep_mixed, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)
tune_by_meansep_del_long <- gather(tune_by_meansep_del, measure, fraction, fdr_all:frac_gnomad, factor_key=TRUE)

pdf('plots/tune_by_mean_sep.by_sv_type.pdf', width = 8, height = 6)
ggplot(tune_by_meansep_dup_long,aes(x=mean_sep_thred, y=fraction, color=measure))+geom_point()+ggtitle("DUP")
ggplot(tune_by_meansep_del_long,aes(x=mean_sep_thred, y=fraction, color=measure))+geom_point()+ggtitle("DEL")
ggplot(tune_by_meansep_mixed_long,aes(x=mean_sep_thred, y=fraction, color=measure))+geom_point()+ggtitle("MIXED")
dev.off()

dup <- filter_callset(var_type = "DUP", MEAN_SEP_thred = 0.47)
del<- filter_callset(var_type = "DEL", MEAN_SEP_thred = 0.47)
mixed<- filter_callset(var_type = "MIXED")
all <- as.data.frame(rbind(dup,del,mixed))
get_fdr(all, min_prob = 2)
