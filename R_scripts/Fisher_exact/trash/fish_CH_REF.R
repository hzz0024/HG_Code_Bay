#######################
#  Fisher's exact test#
#######################
# load reference file with header in it
REF = 'REF_maf0.05_pctind0.7_cv30.mafs'
dat_REF <- read.delim(REF, header = TRUE, sep='\t')
REF_total = dat_REF$nInd*2
REF_refCount = round(dat_REF$knownEM*REF_total)
REF_altCount = REF_total - round(dat_REF$knownEM*REF_total)

# load challenge file with header in it
CH = 'CH_maf0.05_pctind0.7_cv30.mafs'
dat_CH <- read.delim(CH, header = TRUE, sep='\t')
CH_total = dat_CH$nInd*2
CH_refCount = round(dat_CH$knownEM*CH_total)
CH_altCount = CH_total - round(dat_CH$knownEM*CH_total)

outputfile='fish_ch_ref_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(ref_n))
delta_p <- c()
p_values <- c()
alternative <- c()
for(i in idxs){
#for(i in seq(322,322)){
  # print the running process
  s = paste0(i,'/',dim(dat_ref)[1])
  message(s,"\r",appendLF=FALSE)
  # create the matrix for Fisher's exact test
  chr <- dat_ref$chromo[i]
  pos <- dat_ref$position[i]
  c.a1 = REF_refCount[i]
  s.a1 = CH_refCount[i]
  c.a2 = REF_altCount[i]
  s.a2 = CH_altCount[i]
  st = CH_total[i]
  ct = REF_total[i]
  survival = c(s.a1, s.a2)
  control = c(c.a1, c.a2)
  M = as.table(cbind(control, survival))
  rownames(M) = c("ref","alt")
  colnames(M) = c("control", "survival")
  #print(M)
  # Calculate delta_p
  delta_p <- dat_CH$knownEM[i] - dat_REF$knownEM[i]
  #print(delta_p)
  #get T/F values for whether change in allele frequency is positive or negative 
  alternative = delta_p <= 0
  alternative[alternative == TRUE] <- 'greater'
  alternative[alternative == FALSE] <- 'less'
  # perform fisher's exact test
  p = fisher.test(M, alternative=alternative)$p.value
  p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
  p2[p2 > 1] <- 1
  #p_value <- fisher.test(M)$p
  cat(paste0(chr,'\t',pos,'\t',delta_p,'\t',p, '\t',p2,'\n'))
}
sink()

#######################
#  Adjust the p-value #
#######################

library(qvalue)
res_name = 'fish_CH_REF_pvalue.txt'
res <- read.delim(res_name, header = FALSE, sep='\t')
p_value <- res$V5

length(which(p_value<0.05))
hist(p_value)
qobj <- qvalue(p = p_value)
fdr <- qobj$qvalues

length(fdr[fdr<0.2])
outlier <- as.data.frame(cbind(res$V1[fdr<0.2], res$V2[fdr<0.2], p_value[fdr<0.2], fdr[fdr<0.2]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results/Fish_CH_REF_fdr2.txt", quote = FALSE, row.names = FALSE)

length(fdr[fdr<0.1])
outlier <- as.data.frame(cbind(res$V1[fdr<0.1], res$V2[fdr<0.1], p_value[fdr<0.1], fdr[fdr<0.1]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results/Fish_CH_REF_fdr1.txt", quote = FALSE, row.names = FALSE)

length(fdr[fdr<0.05])
outlier <- as.data.frame(cbind(res$V1[fdr<0.05], res$V2[fdr<0.05], p_value[fdr<0.05], fdr[fdr<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
# export the dataframe to txt file
write.table(outlier,"./results/Fish_CH_REF_fdr5.txt", quote = FALSE, row.names = FALSE)

