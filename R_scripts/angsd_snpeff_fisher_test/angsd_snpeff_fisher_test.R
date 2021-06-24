setwd("~/Downloads/physalia_adaptation_course-master/05_day5")

#load library
install.packages("dplyr")
library(dplyr)

#open file
snpEff_db<-read.table("04_snpEff/SNP_annotated_formatted.txt", header=F) 
head(snpEff_db)

#simplify the data
snpEff_db<-snpEff_db[,c(1,2,3,9)] 
head(snpEff_db)

#put informative col names
colnames(snpEff_db)<-c("chr","position","id_snp","category")
head(snpEff_db)

#count how many of each variant category
#with Rbase
all_snps_repartition<-as.matrix(table(snpEff_db$category))
all_snps_repartition

#load data
outlier_temp_rda<-read.table("03_outliers/outlier_temp_rda.txt", header=T)
head(outlier_temp_rda)

#use inner-join function to match the two datatable and keep only rows in common
outlier_annotated<-inner_join(outlier_temp_rda,snpEff_db)
head(outlier_annotated)

#count how many of each variant category
#with Rbase
outlier_repartition<-as.matrix(table(outlier_annotated$category))
outlier_repartition

# we can join both matrix of repartition to have a look
joined_repartition1<- read.table("04_snpEff/joined_repartition1.txt", header=F)
colnames(joined_repartition1)<-c("category","n_all_snps","n_oulier")
joined_repartition1
#joined_repartition1 <- plyr::rbind.fill.matrix(as.matrix(table(snpEff_db$category)), as.matrix(table(outlier_annotated$category)))
#joined_repartition1<-as.data.frame(cbind(as.matrix(table(snpEff_db$category)),as.matrix(table(outlier_annotated$category))))
colnames(joined_repartition1)<-c("n_all_snps","n_oulier")
joined_repartition1$category<-as.character(row.names(joined_repartition1))
joined_repartition1

# And we can use Fisher exact test to ask whether one category of variant is enriched in our outlires
total_snp<-dim(snpEff_db)[1]
total_snp
total_outlier<-dim(outlier_annotated)[1]
total_outlier
# for row i (category i)
i=1
joined_repartition1$category[i]
#build the contigency table
cont_tab<-rbind(c(joined_repartition1$n_oulier[i],total_outlier),
                c(joined_repartition1$n_all_snps[i],total_snp))
cont_tab
cont_tab[ which(is.na(cont_tab))]<-0
fish_res<-fisher.test(cont_tab,alternative = "greater")
fish_res

#You can loop over the categories with a loop and store the results
for (i in 1 : length(joined_repartition1$category))
     {
       cat(joined_repartition1$category[i])
       
       #build the contigency table
       cont_tab<-rbind(c(joined_repartition1$n_oulier[i],total_outlier),
                       c(joined_repartition1$n_all_snps[i],total_snp))
       #cat(cont_tab)
       
       #run Fisher test
       cont_tab[ which(is.na(cont_tab))]<-0
       fish_res<-fisher.test(cont_tab,alternative = "greater")
       print(fish_res$p.value)
}



