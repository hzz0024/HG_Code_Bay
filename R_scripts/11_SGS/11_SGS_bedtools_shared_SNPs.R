
###########################
## test on random shares ##
###########################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS/11_SGS_outlier_share")
#system("/Users/HG/Downloads/bedtools2/bin/bedtools intersect -a challenge_FDR_outlier.list.bed -b HC_NB_FDR_outlier.list.bed -wo > SGS_outlier_shared.txt", intern = TRUE)

bed_format <- function(input, output, window){
  #input = "challenge_FDR_outlier.list.bed"
  outlier_temp_SGS<- read.delim(input, header = FALSE, sep='\t')
  colnames(outlier_temp_SGS) <- c("chr", "pos", "end")
  outlier_temp_SGS$id <- paste0(outlier_temp_SGS$chr, "_", outlier_temp_SGS$pos)
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,5,6,4)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("challenge_FDR_outlier.list.bed", "challenge_FDR_outlier.flank.bed", 10000)
bed_format("HC_NB_FDR_outlier.list.bed", "HC_NB_FDR_outlier.flank.bed", 10000)
system("/Users/HG/Downloads/bedtools2/bin/bedtools intersect -a HC_NB_FDR_outlier.flank.bed -b challenge_FDR_outlier.flank.bed > SGS_outlier_shared_flank.txt", intern = TRUE)



head(daf)
daf$id = paste0(daf$chr,'_',daf$pos)
snp_id = as.data.frame(daf$id)
library(dplyr)
n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  o1 = sample_n(snp_id, 2979)
  o2 = sample_n(snp_id, 33)
  cnt = length(intersect(o1$`daf$id`, o2$`daf$id`))
  boot_cnt[i] = cnt
}
boot_cnt
p = sum(boot_cnt >= 1)/length(boot_cnt)
p
hist(boot_cnt)

# test if overlaps exceed the expection

bed_format <- function(input, output){
  outlier_temp_SGS<- input
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-0
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

daf<-read.table("ps_Del19_challenge.txt", header = F)
colnames(daf) <- c("chr","pos","freqSurv","freqRef", "delta_p", "ps","outlier")
head(daf)
daf$id = paste0(daf$chr,'_',daf$pos)

n_bootstraps=10000
boot_cnt = rep(0, n_bootstraps)
#a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  idx1 = sort(sample.int(dim(daf)[1], 2932))
  idx2 = sort(sample.int(dim(daf)[1], 2507))
  o1 = daf[idx1,]
  o2 = daf[idx2,]
  #o1 = sample_n(daf, 2979)
  #o2 = sample_n(daf, 226)
  bed_format(o1, "test1.bed")
  bed_format(o2, "test2.bed")
  #system("bedtools intersect -a test1.bed -b test2.bed > test.intersect", intern = TRUE)
  #cnt = system("cat test.intersect | wc -l", intern = TRUE)
  system("/Users/HG/Downloads/bedtools2/bin/bedtools intersect -a test1.bed -b test2.bed -wo > test.intersect", intern = TRUE)
  if (file.size("test.intersect") == 0) next
  temp <- read.delim("test.intersect", header = FALSE, sep='\t')
  cnt <- length(unique(c(temp$V4,temp$V8)))
  print(cnt)
  boot_cnt[i] = cnt
}

boot_cnt <- as.numeric(as.character(boot_cnt))
hist(boot_cnt, main="Randomization test for \nshared outliers (n=10000 times)", xlab="Number of SNPs shared", ylab="Frequency")
abline(v = mean(boot_cnt), col = "blue", lwd = 2)
abline(v = 8, col = "red", lwd = 2, lty=2)
text(10, 1500, "observed shared SNPs: 8")
text(10, 1000, "p-value = 0.0546")
# This test confirms whether the normal distribution of the data is violated.
shapiro.test(boot_cnt)
p = sum(boot_cnt >= 8)/length(boot_cnt)
p
t.test(boot_cnt,conf.level=0.95)


###########################
##### find intersects #####
###########################
setwd("~/Documents/Ryan_workplace/DelBay_adult/22_find_intersect")
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/GE_association")
bed_format <- function(input, output){
  outlier_temp_SGS<- read.delim(input, header = TRUE, sep='\t')
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("Del19_FDR_outlier.list", "Del19_FDR_2K.bed")
bed_format("Del20_FDR_outlier.list", "Del20_FDR_2K.bed")
bed_format("NB_HC_FDR_outlier.list", "NB_HC_FDR_2K.bed")
bed_format("SR_HC_FDR_outlier.list", "SR_HC_FDR_2K.bed")

bed_format1 <- function(input, output){
  outlier_temp_SGS<- read.delim(input, header = TRUE, sep='\t')
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$position-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$position+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,17,18,5)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format1("GEA_BF20_outlier.txt", "GEA_BF20_10K.bed")

system("bedtools intersect -a GEA_BF20_10K.bed -b Del19_FDR_10K.bed -wo > GEA_BF_20_Del19_FDR_10K.intersect", intern = TRUE)
cnt = system("cat GEA_BF_20_Del19_FDR_10K.intersect | wc -l", intern = TRUE)
temp <- read.delim("GEA_BF_20_Del19_FDR_10K.intersect", header = FALSE, sep='\t')
cnt <- length(unique(c(temp$V4,temp$V8)))
cnt
