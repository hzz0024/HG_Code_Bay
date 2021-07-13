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
setwd("~/Documents/Ryan_workplace/DelBay_adult/13_env_gen_association/plot_outlier_trend")
setwd("~/Documents/HG/DelBay_adult/13_env_gen_association/salinity2/annotation")

bed_format <- function(input, output){
  #input= "GEA_BF20_outlier.txt"
  outlier_temp_SGS<-read.table(input, header=T)
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
  #outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,18,19,5)]
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,17,18,5)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("GEA_BF20_outlier.txt", "./bed/GEA_BF20_outlier.bed" )
bed_format("GEA_BF10_outlier.txt", "./bed/GEA_BF10_outlier.bed" )
