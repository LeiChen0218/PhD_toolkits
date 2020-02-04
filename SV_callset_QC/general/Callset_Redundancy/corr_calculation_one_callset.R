args <- commandArgs(trailingOnly = T)
if(length(args) < 2)
  stop("input file or output file is missing")

callset_file <- args[1]
if(length(callset_file) < 1)
  stop("input file (callset) is missing")

print("start")
print(Sys.time())

library(reshape2)

callset <- read.table(callset_file, header = T, na.strings = ".")
missing_gt <- apply(callset, 2, function(x) sum(is.na(x))) 
callset <- callset[,missing_gt <= 10]

var_sd <- apply(callset[,-1], 2, function(x) sd(x,na.rm = T))
var_var <- names(var_sd[var_sd != 0])
cn_corr <- cor(callset[var_var],use="pairwise.complete.obs")

cn_corr[is.na(cn_corr)] <- 1

print("finished calculating correlation matrix")
print(Sys.time())

output <- args[2]
write.table(cn_corr, output, col.names = F, row.names = F, quote = F, sep = "\t")
#write.table(cn_corr, 'corr.test',col.names = F, row.names = F, quote = F, sep = "\t")
