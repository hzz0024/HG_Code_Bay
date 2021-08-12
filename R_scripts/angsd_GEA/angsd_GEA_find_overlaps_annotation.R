###############################################################
### Bedtools : find the intersection between SNPs and genes ###
###############################################################

# the code is adopted from the tutorial here https://github.com/clairemerot/physalia_adaptation_course/tree/2021/05_day5

###########################################
### format outlier files into bed files ###
###########################################

#######################
#  find overlaps      #
#######################

# this part is to find the overlaps between Baypass and RDA results (for single SNP)
setwd("~/Documents/HG/DelBay_adult/13_env_gen_association/shared_snp")
outlier_Baypass<-read.table("GEA_BF10_outlier.bed", header=F)
Baypass_outlier_snp <- outlier_Baypass$V4
print(paste0("Baypass BF > 10 has ",length(Baypass_outlier_snp), " outliers"))
outlier_RDA<-read.table("salinity2_env2.txt_candidate_SNP_2_sd.bed", header=F)
RDA_outlier_snp <- outlier_RDA$V4
print(paste0("RDA 2 sd has ",length(RDA_outlier_snp), " outliers"))


# check shared single snp outliers (method1)
Baypass_outlier_snp %in% RDA_outlier_snp
sum(Baypass_outlier_snp %in% RDA_outlier_snp)
# check shared single snp outliers (method2)
length(intersect(Baypass_outlier_snp,RDA_outlier_snp))
# 18

###########################
## test on randomness    ##
###########################
daf<-read.table("ps_Del19_challenge.txt", header = F)
head(daf)
colnames(daf) <- c("chr","pos","freqSurv","freqRef", "delta_p", "ps","outlier")
head(daf)
daf$id = paste0(daf$chr,'_',daf$pos)
snp_id = as.data.frame(daf$id)

library(dplyr)

n_bootstraps=1000 # number of repeats
boot_cnt = rep(0, n_bootstraps)
for (i in 1:n_bootstraps){
  o1 = sample_n(snp_id, 1083, replace=FALSE) # a total of 2032113 snp_id here
  o2 = sample_n(snp_id, 1472, replace=FALSE)
  cnt = length(intersect(o1$`daf$id`, o2$`daf$id`))
  boot_cnt[i] = cnt
}
boot_cnt
median(boot_cnt)
p = sum(boot_cnt >= 18)/length(boot_cnt)
p
hist(boot_cnt)

#######################
#  make bed file      #
#######################
# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("~/Documents/Ryan_workplace/DelBay_adult/13_env_gen_association/plot_outlier_trend")
setwd("~/Documents/HG/DelBay_adult/13_env_gen_association/salinity2/annotation")

bed_format <- function(input, output){
  #input= "GEA_BF20_outlier.txt"
  outlier_temp_SGS<-read.table(input, header=T)
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-1000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$position-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$position+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  #outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,18,19,5)]
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,17,18,5)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("GEA_BF20_outlier.txt", "./bed/GEA_BF20_outlier.bed" )
bed_format("GEA_BF10_outlier.txt", "./bed/GEA_BF10_outlier.bed" )
