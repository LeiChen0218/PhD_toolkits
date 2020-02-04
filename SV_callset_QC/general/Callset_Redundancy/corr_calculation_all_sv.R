#setwd('/Users/leichen/Desktop/Lab/Finmetseq_paper/1-Callset_QC/data')

args <- commandArgs(trailingOnly = T)
if(length(args) < 4)
  stop("input file or output file is missing")

gs_file <- args[1]
if(length(gs_file) < 1)
  stop("input file (gs) is missing")
lumpy_file <- args[2]
if(length(lumpy_file) < 1)
  stop("input file (lumpy) is missing")
cnvnator_file <- args[3]
if(length(cnvnator_file) < 1)
  stop("input file (cnvnator) is missing")

print("start")
print(Sys.time())


library(reshape2)

gs <- read.table(gs_file, header = T, na.strings = ".")
lumpy <- read.table(lumpy_file, header = T, na.strings = ".")
missing_gt <- apply(lumpy, 2, function(x) sum(is.na(x))) 
lumpy <- lumpy[,missing_gt <= 10]
cnvnator <- read.table(cnvnator_file, header = T, na.strings = ".")
#cnvnator <- read.table('CNVnator/highQ_anno/cnvnator.highQ.c21.t.dosage',header = T, na.strings = ".")
#gs <- read.table('GS/highQ_anno/gs.highQ.c21.t.dosage', header = T, na.strings = ".")
#lumpy <- read.table('Lumpy/highQ_anno/lumpy.highQ.c21.t.dosage', header = T, na.strings = ".")

df <- merge(cnvnator,gs,by="ID")
df <- merge(df, lumpy, by="ID")

var_sd <- apply(df[,-1], 2, function(x) sd(x,na.rm = T))
var_var <- names(var_sd[var_sd != 0])
cn_corr <- cor(df[var_var],use="pairwise.complete.obs")

cn_corr[is.na(cn_corr)] <- 1

print("finished calculating correlation matrix")
print(Sys.time())

output <- args[4]
write.table(cn_corr, output, col.names = F, row.names = F, quote = F, sep = "\t")
#write.table(cn_corr, 'corr.test',col.names = F, row.names = F, quote = F, sep = "\t")
