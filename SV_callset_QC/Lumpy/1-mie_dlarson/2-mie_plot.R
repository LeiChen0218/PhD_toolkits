setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/Lumpy/mie')
library(ggplot2)
library(dplyr)

count_per_fam <- read.table('mie_rate_per_trio.txt',header = F)
colnames(count_per_fam) <- c("trio","type","filter","mie_rate")
count_per_fam$type <- as.character(count_per_fam$type)
count_per_fam$type <- factor(count_per_fam$type, levels = c("DEL","MEI","DUP","INV","BND"))
stand_rate <- mean(count_per_fam[which(count_per_fam$type=="DEL"& count_per_fam$filter=="PASS"),]$mie_rate)

ggplot(count_per_fam, aes(x=trio,y=mie_rate))+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  geom_bar(stat="identity")+
  facet_grid(. ~ type)+
  ggtitle('~5k Finn, LUMPY, All SVs')
ggsave('plots/MIE_by_type.png',width = 6,height = 4)

read.table("mie_counts.all_trios.txt", sep="\t", header=T, na.strings=c('.','NA')) -> counts
counts$SVType <- as.character(counts$SVType)
counts$SVType <- as.factor(counts$SVType)

#counts %>% group_by(SVType, MSQBin) %>% summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% mutate(MIERate=MIE/Number) -> counts
ggplot(counts, aes(MSQBin, MIERate, size=Number)) + geom_point() + 
  facet_grid(SVType ~ .)+xlim(0,1000)
ggsave("MSQ_all.png", width=11, height=8.5)

#read.table("mie_counts.all_trios.txt", sep="\t", header=T, na.strings=c('.','NA')) -> counts
#counts %>% group_by(SVType, MSQBin, Filter) %>% summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% mutate(MIERate=MIE/Number) -> counts
#ggplot(counts,  aes(MSQBin, MIERate, size=Number)) + geom_point() + facet_grid(SVType ~ Filter)
#ggsave("MSQ.png", width=11, height=8.5)
#ggplot(counts,  aes(MSQBin, MIERate, size=Number)) + geom_point() + facet_grid(SVType ~ Filter) + xlim(0,1000)
#ggsave("MSQ_zoom.png", width=11, height=8.5)

counts %>% group_by(SVType, MSQBin) %>% 
  summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% 
  arrange(desc(MSQBin)) %>% 
  mutate(cumNumber=cumsum(Number), cumMIE=cumsum(MIE), cumNoMIE=cumsum(NoMIE), cumMIERate=cumMIE/cumNumber) -> x
ggplot(x, aes(MSQBin, cumMIERate)) + geom_point() + 
  facet_grid(SVType ~.) + 
  scale_x_continuous(limits=c(0,1000), minor_breaks = seq(0,1000,25), breaks=seq(0,1000,50)) + 
  ylim(0,0.3) + theme(axis.text.x = element_text(angle=90))+
  ggtitle('Mendelian inheritance Error (MIE) rate  v.s. MSQ')
ggsave("plots/MSQ_cutoff_zoom.png", width=11, height=8.5)

mei <- counts[counts$SVType == "MEI",]
mei %>% group_by(SVType, MSQBin) %>% 
  summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE))%>% 
  arrange(desc(MSQBin)) %>%
  mutate(cumNumber=cumsum(Number), cumMIE=cumsum(MIE), cumNoMIE=cumsum(NoMIE), cumMIERate=cumMIE/cumNumber) -> x

tail(x)
# With PASS variants only
#read.table("mie_counts.all_trios.txt", sep="\t", header=T, na.strings=c('.','NA')) -> counts
#counts <- counts[which(counts$Filter=="PASS"),]

#counts %>% group_by(SVType, MSQBin) %>% summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% arrange(desc(MSQBin)) %>% mutate(cumNumber=cumsum(Number), cumMIE=cumsum(MIE), cumNoMIE=cumsum(NoMIE), cumMIERate=cumMIE/cumNumber) -> x
#ggplot(x, aes(MSQBin, cumMIERate)) + geom_point() + facet_grid(SVType ~ .) + scale_x_continuous(limits=c(0,1000), minor_breaks = seq(0,1000,25), breaks=seq(0,1000,50)) + ylim(0,0.3) + theme(axis.text.x = element_text(angle=90))
#ggsave("MSQ_cutoff_PASS_zoom.png", width=11, height=8.5)

#counts %>% group_by(SVType, MSQBin) %>% summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% arrange(desc(MSQBin)) %>% mutate(cumNumber=cumsum(Number), cumMIE=cumsum(MIE), cumNoMIE=cumsum(NoMIE), cumMIERate=cumMIE/cumNumber) -> x
#ggplot(x, aes(MSQBin, cumMIERate)) + geom_point() + facet_grid(SVType ~ .) + theme(axis.text.x = element_text(angle=90))
#ggsave("MSQ_cutoff_LOW_fullscale_zoom.png", width=11, height=8.5)


# With LOW variants only
#read.table("mie_counts.all_trios.txt", sep="\t", header=T, na.strings=c('.','NA')) -> counts
#counts <- counts[which(counts$Filter=="LOW"),]

counts %>% group_by(SVType, MSQBin) %>% summarize(Number=sum(Number), MIE=sum(MIE), NoMIE=sum(NoMIE)) %>% arrange(desc(MSQBin)) %>% mutate(cumNumber=cumsum(Number), cumMIE=cumsum(MIE), cumNoMIE=cumsum(NoMIE), cumMIERate=cumMIE/cumNumber) -> x
ggplot(x, aes(MSQBin, cumMIERate)) + geom_point() + facet_grid(SVType ~ .) + scale_x_continuous(limits=c(0,1000), minor_breaks = seq(0,1000,25), breaks=seq(0,1000,50)) + ylim(0,0.3) + theme(axis.text.x = element_text(angle=90))
ggsave("MSQ_cutoff_LOW_zoom.png", width=11, height=8.5)