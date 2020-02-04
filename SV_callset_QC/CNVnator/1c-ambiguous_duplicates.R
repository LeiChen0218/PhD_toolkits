setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/CNVnator/IRS_tuning')
library(ggplot2)
library(tidyr)

common <- read.table('common.tuning.table', header = T)
rare <- read.table('rare.tuning.table', header = T)
dups <- read.table('duplicates.list')
colnames(dups) <- c("ID","common_info","rare_info")

dups <- merge(dups, common, by = "ID")
dups <- merge(dups, rare, by = "ID")
head(dups)

ggplot(dups, aes(x=DIP_P, y=NNONREF))+geom_point()
ggplot(dups, aes(x=DIP_P))+geom_histogram(bins=50)+
  ggtitle("dip p distribution on ~600 ambiguous duplicates")
ggsave('plots/dip_p_dist.amb_dup.png',width = 6, height = 4)

write.table(dups[dups$DIP_P > 0.05,], "remove_from_common.txt",
            quote = F, row.names = F, sep = "\t")

write.table(dups[dups$DIP_P <= 0.05,], "remove_from_rare.txt",
            quote = F, row.names = F, sep = "\t")
