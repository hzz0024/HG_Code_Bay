setwd("~/Downloads/physalia_adaptation_course-master/05_day5")

#load file
outlier_temp_rda<-read.table("03_outliers/outlier_temp_rda.txt", header=T)
#have a quick look
head(outlier_temp_rda)
#what's its dimension?
dim(outlier_temp_rda)

#which size around the SNP
window<-10000
#add a vector with start position
outlier_temp_rda$start<-outlier_temp_rda$position-(window/2)
#oups it can't be ngative! replace negative by 0
outlier_temp_rda$start[outlier_temp_rda$start<0]<-0 
#add a vector with stop position
outlier_temp_rda$stop<-outlier_temp_rda$position+(window/2)
#have a look
head(outlier_temp_rda)
#which columns shoud we keep?
outlier_temp_rda_bed<-outlier_temp_rda[,c(2,5,6,1)]
#save your file
write.table(outlier_temp_rda_bed, "05_bed/outlier_temp_rda.bed", row.names=F, sep="\t", quote=F,col.names=F)

all_snps<-read.table("03_outliers/SNP_pos.txt", header=T)
all_snps$start<-all_snps$position-(window/2)
all_snps$start[all_snps$start<0]<-0 
all_snps$stop<-all_snps$position+(window/2)
head(all_snps)
all_snps_bed<-all_snps[,c(2,4,5,1)]
write.table(all_snps_bed, "05_bed/all_snps.bed", row.names=F, sep="\t", quote=F,col.names=F)




