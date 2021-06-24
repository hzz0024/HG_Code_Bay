###############################################################
### Bedtools : find the intersection between SNPs and genes ###
###############################################################

# the code is adopted from the tutorial here https://github.com/clairemerot/physalia_adaptation_course/tree/2021/05_day5

###########################################
### format outlier files into bed files ###
###########################################


#######################
#  make bed file      #
#######################
# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("~/Documents/Ryan_workplace/DelBay_adult/11_SGS")

bed_format <- function(input, output){
  outlier_temp_SGS<-read.table(input, header=T)
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

bed_format("Del19_FDR_outlier.list", "./bed/Del19_FDR_10K.bed" )
bed_format("Del20_FDR_outlier.list", "./bed/Del20_FDR_10K.bed" )
bed_format("NB_HC_FDR_outlier.list", "./bed/NB_HC_FDR_10K.bed" )
bed_format("SR_HC_FDR_outlier.list", "./bed/SR_HC_FDR_10K.bed" )

# do the same thing for whole SNPs

outlier_temp_SGS<-read.table("ps_Del19_challenge.txt", header=F)
#what's its dimension?
dim(outlier_temp_SGS)
head(outlier_temp_SGS)
#which size around the SNP
window<-10000
#add a vector with start position
outlier_temp_SGS$start<-outlier_temp_SGS$V2-(window/2)
#oups it can't be ngative! replace negative by 0
outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
#add a vector with stop position
outlier_temp_SGS$stop<-outlier_temp_SGS$V2+(window/2)
outlier_temp_SGS$id_snp <- paste0(outlier_temp_SGS$V1,'_',outlier_temp_SGS$V2)
#have a look
head(outlier_temp_SGS)
#which columns shoud we keep?
outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,8,9,10)]
colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
head(outlier_temp_SGS_bed)
#save your file
write.table(outlier_temp_SGS_bed, "./bed/all_snps.bed", row.names=F, sep="\t", quote=F,col.names=F)

###########################################
#  make bed file for Fisher outliers      #
###########################################

# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("~/Documents/Ryan_workplace/DelBay_adult/08_fish_exact/annotation_z_outlier")

bed_format <- function(input, output){
  outlier_temp_SGS<-read.table(input, header=F)
  #outlier_temp_SGS<-read.table("REF19_CHR19_NB_HC_out_0.05_fish.txt", header=F)
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  #window<-2000
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$V2-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$V2+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,4,5,3)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "FDR")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("REF19_CHR19_NB_HC_out_0.05_fish.txt", "./bed/REF19_CHR19_NB_HC_out_0.05_2K.bed" )
bed_format("REF19_CHR19_NB_HC_out_0.05_fish.txt", "./bed/REF19_CHR19_NB_HC_out_0.05_10K.bed" )
bed_format("REF19_CHR19_SR_HC_out_0.05_fish.txt", "./bed/REF19_CHR19_SR_HC_out_0.05_2K.bed" )
bed_format("REF19_CHR19_SR_HC_out_0.05_fish.txt", "./bed/REF19_CHR19_SR_HC_out_0.05_10K.bed" )


setwd("~/Documents/Ryan_workplace/DelBay_adult/08_fish_exact/annotation_z_outlier/bed/")
system("bedtools intersect -a REF19_CHR19_NB_HC_out_0.05_2K.bed -b REF19_CHR19_SR_HC_out_0.05_2K.bed -wo > REF19_CHR19_SR_HC_NB_HC_FDR_2K.intersect", intern = TRUE)
system("bedtools intersect -a REF19_CHR19_NB_HC_out_0.05_10K.bed -b REF19_CHR19_SR_HC_out_0.05_10K.bed -wo > REF19_CHR19_SR_HC_NB_HC_FDR_10K.intersect", intern = TRUE)
