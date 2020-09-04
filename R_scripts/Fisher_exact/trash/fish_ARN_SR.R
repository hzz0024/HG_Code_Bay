#######################
#  Fisher's exact test#
#######################
ARN = 'ARN_maf0.05_pctind0.7_cv30.mafs'
dat_ARN <- read.delim(ARN, header = TRUE, sep='\t')
ARN_n = dat_ARN$nInd
ARN_k = round(dat_ARN$knownEM*dat_ARN$nInd*2)

SR = 'SR_maf0.05_pctind0.7_cv30.mafs'
dat_SR <- read.delim(SR, header = TRUE, sep='\t')
SR_n <- dat_SR$nInd
SR_k <- round(dat_SR$knownEM*dat_SR$nInd*2)

outputfile='fish_ARN_SR_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(ARN_n))
p_values <- c()
for(i in idxs){
#for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_ARN)[1])
  message(s,"\r",appendLF=FALSE)

  chr <- dat_ARN$chromo[i]
  pos <- dat_ARN$position[i]
  mat <- matrix(c(SR_k[i],(2*SR_n[i])-SR_k[i],ARN_k[i],(2*ARN_n[i])-ARN_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()