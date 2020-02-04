setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data/general/genome_density')
library('ggplot2')
library(tidyr)

coord <- read.table('b38.1mb.windows.bed')
colnames(coord) <- c("chr","start","end")
coord$genome <- as.numeric(rownames(coord))

# gs
gs <- read.table('gs.1mb.density.bed')
colnames(gs) <- c("chr","start","end","gs")
coord <- merge(coord, gs, by = c("chr","start","end"), all.x = T)

# lumpy
lumpy <- read.table('lumpy.1mb.density.bed')
colnames(lumpy) <- c("chr","start","end","lumpy")
coord <- merge(coord, lumpy, by = c("chr","start","end"), all.x = T)

# cnvnator
cnvnator <- read.table('cnvnator.1mb.density.bed')
colnames(cnvnator) <- c("chr","start","end","cnvnator")
coord <- merge(coord, cnvnator, by = c("chr","start","end"), all.x = T)

coord$color <- "odd"
coord[coord$chr %in% c("chr2","chr4","chr6","chr8","chr10","chr12","chr14","chr16","chr18","chr20","chr22"),]$color <- "even"

# reformate
coord[is.na(coord)] <- 0
df <- gather(coord, callset, count, gs:cnvnator)

ggplot(df, aes(x=genome, y=count, fill=color))+geom_bar(stat="identity")+
  facet_grid(callset~., scales = "free_y")+ylab('SV per-Mb')+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab('genome')+
  theme(legend.position="none")+ggtitle('high-confident SV density on each chromosome ')
ggsave('plots/variant_density.per_chr.png', width = 10, height = 6)

ggplot(df, aes(x=count,fill=callset))+geom_histogram(bins = 300)+
  facet_grid(callset~., scales = "free_y")+xlab('variant count in 1MB windows')
ggsave('plots/variant_density.hist.png')

