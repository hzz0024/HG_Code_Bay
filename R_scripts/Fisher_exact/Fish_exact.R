#######################
#  Fisher's exact test#
#######################
# load reference file with header in it
ref = 'REF_maf0.05_pctind0.7_cv30.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
ref_n = dat_ref$nInd
ref_k = round(dat_ref$knownEM*dat_ref$nInd*2)

# load challenge file with header in it
chr = 'CH_maf0.05_pctind0.7_cv30.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
ch_n <- dat_ch$nInd
ch_k <- round(dat_ch$knownEM*dat_ch$nInd*2)

outputfile='RES.txt'
sink(outputfile)
idxs <- seq(1, length(ref_n))
p_values <- c()
for(i in idxs){
#for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_ref)[1])
  message(s,"\r",appendLF=FALSE)
  chr <- dat_ref$chromo[i]
  pos <- dat_ref$position[i]
  mat <- matrix(c(ref_k[i],(2*ref_n[i])-ref_k[i],ch_k[i],(2*ch_n[i])-ch_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()

#######################
#  Adjust the p-value #
#######################
# # load SNP dataset
library(qvalue)
res_name = 'RES.txt'
res <- read.delim(res_name, header = FALSE, sep='\t')
p_value <- res$V3

length(which(p_value<0.01))
hist(p_value)
p_value[p_value > 1] = 1
qobj <- qvalue(p = p_value)
lfdr <- qobj$lfdr
length(lfdr[lfdr<0.05])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.05], res$V2[lfdr<0.05], p_value[lfdr<0.05], lfdr[lfdr<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to csv file
write.csv(outlier,"fisher_exact_outlier.csv", row.names = FALSE)
