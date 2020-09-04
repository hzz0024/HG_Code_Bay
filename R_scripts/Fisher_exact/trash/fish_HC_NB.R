#######################
#  Fisher's exact test#
#######################
HC = 'HC_maf0.05_pctind0.7_cv30.mafs'
dat_HC <- read.delim(HC, header = TRUE, sep='\t')
HC_n = dat_HC$nInd
HC_k = round(dat_HC$knownEM*dat_HC$nInd*2)

NB = 'NB_maf0.05_pctind0.7_cv30.mafs'
dat_NB <- read.delim(NB, header = TRUE, sep='\t')
NB_n <- dat_NB$nInd
NB_k <- round(dat_NB$knownEM*dat_NB$nInd*2)

outputfile='fish_NB_HC_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(HC_n))
p_values <- c()
for(i in idxs){
#for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_HC)[1])
  message(s,"\r",appendLF=FALSE)

  chr <- dat_HC$chromo[i]
  pos <- dat_HC$position[i]
  mat <- matrix(c(NB_k[i],(2*NB_n[i])-NB_k[i],HC_k[i],(2*HC_n[i])-HC_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()
