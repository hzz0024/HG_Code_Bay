name = 'p_value_local_100k.txt'
pvalue <- read.delim(name, header = FALSE, sep='\t')
cnt <- length(pvalue$V4[pvalue$V4 == 1])
length(pvalue$V3[pvalue$V3 == 0])
hist(pvalue$V3)

name = 'p_value_local_1m.txt'
pvalue <- read.delim(name, header = FALSE, sep='\t')
cnt <- length(pvalue$V4[pvalue$V4 == 1])
length(pvalue$V3[pvalue$V3 == 0])
hist(pvalue$V3)


#######################
#  Adjust the p-value #
#######################
# load SNP dataset
library(qvalue)
name = 'p_value_local_1m.txt'
res <- read.delim(name, header = FALSE, sep='\t')
p_value <- res$V3
length(which(p_value<0.05))
hist(p_value)
qobj <- qvalue(p = p_value,pi0=1)
lfdr <- qobj$lfdr
length(lfdr[lfdr<0.05])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.5], res$V2[lfdr<0.5], p_value[lfdr<0.5], lfdr[lfdr<0.5]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to csv file
write.table(outlier,"SGS_local_1M_outlier.txt", row.names = FALSE, sep="\t", quote = FALSE)


# load SNP dataset
library(qvalue)
name = 'p_value_local_100k.txt'
res <- read.delim(name, header = FALSE, sep='\t')
p_value <- res$V3
length(which(p_value<0.05))
hist(p_value)
qobj <- qvalue(p = p_value,pi0=1)
lfdr <- qobj$lfdr
length(lfdr[lfdr<0.05])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.5], res$V2[lfdr<0.5], p_value[lfdr<0.5], lfdr[lfdr<0.5]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to csv file
write.table(outlier,"SGS_local_100k_outlier.txt", row.names = FALSE, sep="\t", quote = FALSE)
