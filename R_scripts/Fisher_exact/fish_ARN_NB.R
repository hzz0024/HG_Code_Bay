#######################
#  Fisher's exact test#
#######################
ARN = 'ARN_maf0.05_pctind0.7_cv30.mafs'
dat_ARN <- read.delim(ARN, header = TRUE, sep='\t')
ARN_n = dat_ARN$nInd
ARN_k = round(dat_ARN$knownEM*dat_ARN$nInd*2)

NB = 'NB_maf0.05_pctind0.7_cv30.mafs'
dat_NB <- read.delim(NB, header = TRUE, sep='\t')
NB_n <- dat_NB$nInd
NB_k <- round(dat_NB$knownEM*dat_NB$nInd*2)

outputfile='fish_ARN_NB_pvalue.txt'
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
  mat <- matrix(c(NB_k[i],(2*NB_n[i])-NB_k[i],ARN_k[i],(2*ARN_n[i])-ARN_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()