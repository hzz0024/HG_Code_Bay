#######################
#  Fisher's exact test#
#######################
COH = 'COH_maf0.05_pctind0.7_cv30.mafs'
dat_COH <- read.delim(COH, header = TRUE, sep='\t')
COH_n = dat_COH$nInd
COH_k = round(dat_COH$knownEM*dat_COH$nInd*2)

SR = 'SR_maf0.05_pctind0.7_cv30.mafs'
dat_SR <- read.delim(SR, header = TRUE, sep='\t')
SR_n <- dat_SR$nInd
SR_k <- round(dat_SR$knownEM*dat_SR$nInd*2)

outputfile='fish_COH_SR_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(COH_n))
p_values <- c()
for(i in idxs){
#for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_COH)[1])
  message(s,"\r",appendLF=FALSE)

  chr <- dat_COH$chromo[i]
  pos <- dat_COH$position[i]
  mat <- matrix(c(SR_k[i],(2*SR_n[i])-SR_k[i],COH_k[i],(2*COH_n[i])-COH_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()