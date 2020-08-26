library(qvalue)
res_name = 'fish_CH_REF_pvalue.txt'
res <- read.delim(res_name, header = FALSE, sep='\t')
p_value <- res$V3

length(which(p_value<0.05))
hist(p_value)
p_value[p_value > 1] = 1
qobj <- qvalue(p = p_value)
lfdr <- qobj$lfdr

length(lfdr[lfdr<0.2])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.2], res$V2[lfdr<0.2], p_value[lfdr<0.2], lfdr[lfdr<0.2]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_CH_REF_fdr2.txt", quote = FALSE, row.names = FALSE)

length(lfdr[lfdr<0.1])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.1], res$V2[lfdr<0.1], p_value[lfdr<0.1], lfdr[lfdr<0.1]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_CH_REF_fdr1.txt", quote = FALSE, row.names = FALSE)

length(lfdr[lfdr<0.05])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.05], res$V2[lfdr<0.05], p_value[lfdr<0.05], lfdr[lfdr<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_CH_REF_fdr5.txt", quote = FALSE, row.names = FALSE)

library(qvalue)
res_name = 'fish_HC_NB_pvalue.txt'
res <- read.delim(res_name, header = FALSE, sep='\t')
p_value <- res$V3

length(which(p_value<0.05))
hist(p_value)
p_value[p_value > 1] = 1
qobj <- qvalue(p = p_value)
lfdr <- qobj$lfdr

length(lfdr[lfdr<0.2])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.2], res$V2[lfdr<0.2], p_value[lfdr<0.2], lfdr[lfdr<0.2]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_HC_NB_fdr2.txt", quote = FALSE, row.names = FALSE)

length(lfdr[lfdr<0.1])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.1], res$V2[lfdr<0.1], p_value[lfdr<0.1], lfdr[lfdr<0.1]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_HC_NB_fdr1.txt", quote = FALSE, row.names = FALSE)

length(lfdr[lfdr<0.05])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.05], res$V2[lfdr<0.05], p_value[lfdr<0.05], lfdr[lfdr<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_HC_NB_fdr5.txt", quote = FALSE, row.names = FALSE)

# check common shared outliers

CH_REF_name = 'Fish_CH_REF_fdr2.txt'
CH_REF <- read.delim(CH_REF_name, header = TRUE, sep=' ')
CH_REF_list = paste0(CH_REF$chr,'_',CH_REF$pos)
length(CH_REF_list)

HC_NB_name = 'Fish_HC_NB_fdr2.txt'
HC_NB <- read.delim(HC_NB_name, header = TRUE, sep=' ')
HC_NB_list = paste0(HC_NB$chr,'_',HC_NB$pos)
length(HC_NB_list)

CH_REF_name = 'Fish_CH_REF_fdr1.txt'
CH_REF <- read.delim(CH_REF_name, header = TRUE, sep=' ')
CH_REF_list = paste0(CH_REF$chr,'_',CH_REF$pos)
length(CH_REF_list)

HC_ARN_name = 'Fish_HC_ARN_fdr1.txt'
HC_ARN <- read.delim(HC_ARN_name, header = TRUE, sep=' ')
HC_ARN_list = paste0(HC_ARN$chr,'_',HC_ARN$pos)
length(HC_ARN_list)

HC_COH_name = 'Fish_HC_COH_fdr1.txt'
HC_COH <- read.delim(HC_COH_name, header = TRUE, sep=' ')
HC_COH_list = paste0(HC_COH$chr,'_',HC_COH$pos)
length(HC_COH_list)

HC_SR_name = 'Fish_HC_SR_fdr1.txt'
HC_SR <- read.delim(HC_SR_name, header = TRUE, sep=' ')
HC_SR_list = paste0(HC_SR$chr,'_',HC_SR$pos)
length(HC_SR_list)

HC_NB_name = 'Fish_HC_NB_fdr1.txt'
HC_NB <- read.delim(HC_NB_name, header = TRUE, sep=' ')
HC_NB_list = paste0(HC_NB$chr,'_',HC_NB$pos)
length(HC_NB_list)

ARN_COH_name = 'Fish_ARN_COH_fdr1.txt'
ARN_COH <- read.delim(ARN_COH_name, header = TRUE, sep=' ')
ARN_COH_list = paste0(ARN_COH$chr,'_',ARN_COH$pos)
length(ARN_COH_list)

ARN_SR_name = 'Fish_ARN_SR_fdr1.txt'
ARN_SR <- read.delim(ARN_SR_name, header = TRUE, sep=' ')
ARN_SR_list = paste0(ARN_SR$chr,'_',ARN_SR$pos)
length(ARN_SR_list)

ARN_NB_name = 'Fish_ARN_NB_fdr1.txt'
ARN_NB <- read.delim(ARN_NB_name, header = TRUE, sep=' ')
ARN_NB_list = paste0(ARN_NB$chr,'_',ARN_NB$pos)
length(ARN_NB_list)

COH_SR_name = 'Fish_COH_SR_fdr1.txt'
COH_SR <- read.delim(COH_SR_name, header = TRUE, sep=' ')
COH_SR_list = paste0(COH_SR$chr,'_',COH_SR$pos)
length(COH_SR_list)

COH_NB_name = 'Fish_COH_NB_fdr1.txt'
COH_NB <- read.delim(COH_NB_name, header = TRUE, sep=' ')
COH_NB_list = paste0(COH_NB$chr,'_',COH_NB$pos)
length(COH_NB_list)

SR_NB_name = 'Fish_SR_NB_fdr1.txt'
SR_NB <- read.delim(SR_NB_name, header = TRUE, sep=' ')
SR_NB_list = paste0(SR_NB$chr,'_',SR_NB$pos)
length(SR_NB_list)

intersect(HC_ARN_list, HC_COH_list)
intersect(HC_ARN_list, HC_SR_list)
intersect(HC_ARN_list, HC_NB_list)
intersect(ARN_COH_list, ARN_SR_list)
intersect(ARN_COH_list, ARN_NB_list)
intersect(COH_SR_list, COH_NB_list)

a <- intersect(HC_ARN_list, ARN_COH_list)
b <- intersect(ARN_COH_list, COH_SR_list)
c <- intersect(COH_SR_list, SR_NB_list)
intersect(intersect(a, b), c)

intersect(CH_REF_list, HC_COH_list)
intersect(CH_REF_list, HC_ARN_list)
intersect(CH_REF_list, HC_SR_list)
intersect(CH_REF_list, HC_NB_list)
intersect(CH_REF_list, ARN_COH_list)
intersect(CH_REF_list, ARN_SR_list)
intersect(CH_REF_list, ARN_NB_list)
intersect(CH_REF_list, COH_SR_list)
intersect(CH_REF_list, COH_NB_list)
intersect(CH_REF_list, SR_NB_list)


