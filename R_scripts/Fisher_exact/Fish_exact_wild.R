#######################
#  Fisher's exact test#
#######################
# load reference file with header in it
NB = 'NB_maf0.05_pctind0.7_cv30.mafs'
dat_NB <- read.delim(NB, header = TRUE, sep='\t')
NB_n = dat_NB$nInd
NB_k = round(dat_NB$knownEM*dat_NB$nInd*2)

# load challenge file with header in it
HC = 'HC_maf0.05_pctind0.7_cv30.mafs'
dat_HC <- read.delim(HC, header = TRUE, sep='\t')
HC_n <- dat_HC$nInd
HC_k <- round(dat_HC$knownEM*dat_HC$nInd*2)

outputfile='RES_wild.txt'
sink(outputfile)
idxs <- seq(1, length(NB_n))
p_values <- c()
for(i in idxs){
  #for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_NB)[1])
  message(s,"\r",appendLF=FALSE)
  chr <- dat_NB$chromo[i]
  pos <- dat_NB$position[i]
  mat <- matrix(c(HC_k[i],(2*HC_n[i])-HC_k[i],NB_k[i],(2*NB_n[i])-NB_k[i]),nrow=2)
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
res_name = 'RES_wild.txt'
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
# export the dataframe to txt file
write.table(outlier,"fisher_exact_outlier_wild.txt", quote = FALSE, row.names = FALSE)
