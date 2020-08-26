#######################
#  Fisher's exact test#
#######################
SR = 'SR_maf0.05_pctind0.7_cv30.mafs'
dat_SR <- read.delim(SR, header = TRUE, sep='\t')
SR_n = dat_SR$nInd
SR_k = round(dat_SR$knownEM*dat_SR$nInd*2)

NB = 'NB_maf0.05_pctind0.7_cv30.mafs'
dat_NB <- read.delim(NB, header = TRUE, sep='\t')
NB_n <- dat_NB$nInd
NB_k <- round(dat_NB$knownEM*dat_NB$nInd*2)

outputfile='fish_SR_NB_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(SR_n))
p_values <- c()
for(i in idxs){
#for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_SR)[1])
  message(s,"\r",appendLF=FALSE)

  chr <- dat_SR$chromo[i]
  pos <- dat_SR$position[i]
  mat <- matrix(c(NB_k[i],(2*NB_n[i])-NB_k[i],SR_k[i],(2*SR_n[i])-SR_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()