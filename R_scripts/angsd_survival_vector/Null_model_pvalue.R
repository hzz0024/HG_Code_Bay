#######################
#  Produce p-value    #
#######################
# load the allele frequency inputs with ~1.7 million SNPs
ref_name = 'CHR_maf0.05_pctind0.7_cv30.mafs'  
ref_dat = read.delim(ref_name, header = TRUE, sep='\t')
ref_af = ref_dat$knownEM # this is the allele frquency values before selection (i.e. p0 in the slide)
ref_chr = ref_dat$chromo
ref_pos = ref_dat$position

ch_name = 'CH_maf0.05_pctind0.7_cv30.mafs'  
ch_dat = read.delim(ch_name, header = TRUE, sep='\t')
ch_af = ch_dat$knownEM # this is the allele frquency values after selection (i.e. p1 in the slide)

obs_deltaps = ch_af - ref_af # this is observed delta_p

# load the delta_p values from 800 simulations
data = data.frame(matrix(vector(), dim(ref_dat)[1], 800))
for(idx in seq(1,800)){
  message(idx)
  file = paste0('out',idx,'/','test_lrt_obs')
  dat = read.delim(file, header = FALSE, sep=' ')
  MAF = dat$V1 
  data[idx] = MAF - ref_af
}
# Calculate the p_values
outputfile='p_values.txt'
sink(outputfile)
cat(paste0('chr\tpos\tpvalue\n'))
ps = c()
#for(snp in seq(1,500)){   # test for the first 500 SNPs
for(snp in seq(1,dim(ref_dat)[1])){
  s = paste0(snp,'/',dim(ref_dat)[1])
  message(s,"\r",appendLF=FALSE)
  
  obs_deltap = obs_deltaps[snp]
  cat(ref_dat$chromo[snp])
  cat('\t')
  cat(ref_dat$position[snp])
  cat('\t')
  delta_ps = as.numeric(data[snp,])
  # function for p-value calculation, the +1 here is applied to avoid producing p_value == 0
  if(obs_deltap>0){
    p_value <- (length(delta_ps[delta_ps>=obs_deltap])+1)/(length(delta_ps)+1)
  }
  if(obs_deltap<0){
    p_value <- (length(delta_ps[delta_ps<=obs_deltap])+1)/(length(delta_ps)+1)
  }
  if(obs_deltap==0){
    p_value <- 1
  }
  ps = c(ps, p_value)
  cat(p_value)
  cat('\n')
}
sink()

#######################
#  Adjust the p-value #
#######################
# load SNP dataset
library(qvalue)
res_name = 'p_values.txt'
res <- read.delim(res_name, header = TRUE, sep='\t')
p_value <- res$pvalue
# check how many SNPs with p_value < 0.05
length(which(p_value<0.05))
# plot the p-value histogram
hist(p_value)
# start p_value adjustment
qobj <- qvalue(p = p_value)
lfdr <- qobj$lfdr
# check how many SNPs with FDR < 0.2
length(lfdr[lfdr<0.2])
outlier <- as.data.frame(cbind(res$chr[lfdr<0.2], res$pos[lfdr<0.2], p_value[lfdr<0.2], lfdr[lfdr<0.2]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"survival_vector_outlier.txt", row.names = FALSE, sep="\t", quote = FALSE)
