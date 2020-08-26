#######################
#  Fisher's exact test#
#######################
# load reference file with header in it
REF = 'REF_maf0.05_pctind0.7_cv30.mafs'
dat_REF <- read.delim(REF, header = TRUE, sep='\t')
REF_n = dat_REF$nInd
REF_k = round(dat_REF$knownEM*dat_REF$nInd*2)

# load challenge file with header in it
CH = 'CH_maf0.05_pctind0.7_cv30.mafs'
dat_CH <- read.delim(CH, header = TRUE, sep='\t')
CH_n <- dat_CH$nInd
CH_k <- round(dat_CH$knownEM*dat_CH$nInd*2)

outputfile='fish_CH_REF_pvalue.txt'
sink(outputfile)
idxs <- seq(1, length(REF_n))
p_values <- c()
#for(i in idxs){
for(i in seq(1,1000)){
  # print the running process
  s = paste0(i,'/',dim(dat_REF)[1])
  message(s,"\r",appendLF=FALSE)
  
  chr <- dat_REF$chromo[i]
  pos <- dat_REF$position[i]
  p0 <- dat_REF$knownEM[i]
  p1 <- dat_CH$knownEM[i]
  dp <- p1 - p0
  mat <- matrix(c(CH_k[i],(2*CH_n[i])-CH_k[i],REF_k[i],(2*REF_n[i])-REF_k[i]),nrow=2)
  p_value <- fisher.test(mat)$p
  p_values <- c(p_values, p_value)
  cat(paste0(chr,'\t',pos,'\t',p0,'\t',p1,'\t',dp,'\t', p_value,'\n'))
}
sink()

#######################
#  Adjust the p-value #
#######################

library(qvalue)
res_name = 'fish_CH_REF_pvalue.txt'
res <- read.delim(res_name, header = FALSE, sep='\t')
p_value <- res$V3

length(which(p_value<0.05))
hist(p_value)
p_value[p_value > 1] = 1
qobj <- qvalue(p = p_value)
lfdr <- qobj$lfdr
length(lfdr[lfdr<0.1])
outlier <- as.data.frame(cbind(res$V1[lfdr<0.1], res$V2[lfdr<0.1], p_value[lfdr<0.1], lfdr[lfdr<0.1]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_CH_REF_fdr1.txt", quote = FALSE, row.names = FALSE)

outlier <- as.data.frame(cbind(res$V1[lfdr<0.05], res$V2[lfdr<0.05], p_value[lfdr<0.05], lfdr[lfdr<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"./results_fdr/Fish_CH_REF_fdr5.txt", quote = FALSE, row.names = FALSE)

##########################
#  Reveal the relation-  #
#  ship between deltap   #
#and start p, first plot #
##########################

##################### reveal the relationship between deltap and start p, first plot #####################
setwd("/Volumes/cornell/Fisher_exact/results_fdr")
# Delta p plot
library(ggplot2)
# load the neutral datasets
# files <- list.files('.', pattern = "*.txt")
# load the dataset
CH_file = 'CH_maf0.05_pctind0.7_cv30.mafs.extracted'
REF_file = 'REF_maf0.05_pctind0.7_cv30.mafs.extracted'
# obtain the deltap from obs dataset
CH = read.delim(CH_file, header = FALSE, sep = "\t", dec = ".")
REF = read.delim(REF_file, header = FALSE, sep = "\t", dec = ".")
p0 = REF$V6
p1 = CH$V6
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = 96 # change this number according to the outlier file
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNP", limits=c(0, 100)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = "Allele frequency changes for 96 potential outliers identified from Fisher's exact test",
       subtitle = "REF allele = black, CH allele = red, actual deltap = yellow")

##################### reveal the relationship between deltap and start p, second plot #####################
library(ggplot2)
# load the dataset
CH_file = 'CH_maf0.05_pctind0.7_cv30.mafs.extracted'
REF_file = 'REF_maf0.05_pctind0.7_cv30.mafs.extracted'
#deltap_file = 'obs_deltap.output'
# obtain the deltap from obs dataset
CH = read.delim(CH_file, header = FALSE, sep = "\t", dec = ".")
REF = read.delim(REF_file, header = FALSE, sep = "\t", dec = ".")
deltap = CH$V6 - REF$V6
p0 = REF$V6
p1 = CH$V6
DATA = data.frame(p=p0, delta_p=deltap)
DATA = DATA[order(DATA$p),]
num_snp = 96
sp <- ggplot(DATA, aes(x=p, y=delta_p)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p", limits=c(min(DATA$p), max(DATA$p))) +
  scale_y_continuous(name="Deltap (absolute values)", limits=c(min(DATA$delta_p), max(DATA$delta_p))) +
  labs(title = "Deltap against reference allele p for Fisher's exact result (96 outliers after FDR correction)",
       subtitle = "note deltap ranges from 0-0.5")