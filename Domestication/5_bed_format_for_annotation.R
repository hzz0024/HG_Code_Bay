#######################
#  make bed file      #
#######################
# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("~/Documents/HG/Domestication/Annotation/PCAdapt_outliers")

bed_format <- function(input, output){
  outlier_temp<-read.table(input, header=F)
  #outlier_temp<-read.table("Ind514_best_practice_outlier_PC1-2_10K.txt", header=F)
  # replace chromosome if it is numerical
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_temp$V1[outlier_temp$V1==i] = chr_str_list[i]
  outlier_temp[with(outlier_temp, order(V1, V2)),]
  #what's its dimension?
  dim(outlier_temp)
  outlier_temp$V2 <- as.numeric(outlier_temp$V2)
  #which size around the SNP
  window<-0
  #add a vector with start position
  outlier_temp$start<-outlier_temp$V2-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp$start[outlier_temp$start<0]<-0 
  #add a vector with stop position
  outlier_temp$stop<-outlier_temp$V2+(window/2)
  #have a look
  head(outlier_temp)
  #which columns shoud we keep?
  outlier_temp_bed<-outlier_temp[,c(1,5,6,4)]
  colnames(outlier_temp_bed) <- c("chromosome", "start", "stop", "neg_log10")
  head(outlier_temp_bed)
  #save your file
  write.table(outlier_temp_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("Ind514_best_practice_outlier_PC1-2_10K.txt", "./Ind514_best_practice_outlier_PC1-2_10K.bed")
bed_format("Ind514_best_practice_outlier_no_chr156invers_PC1-2_10K.txt", "./Ind514_best_practice_outlier_no_chr156invers_PC1-2_10K.bed")



