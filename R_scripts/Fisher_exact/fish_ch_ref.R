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

outputfile='fish_ch_ref_pvalue.txt'
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
  mat <- matrix(c(ch_k[i],(2*ch_n[i])-ch_k[i],ref_k[i],(2*ref_n[i])-ref_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p_value,'\n'))
}
sink()
