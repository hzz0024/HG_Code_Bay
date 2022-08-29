setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/21_ngsparalog/validation")
library(ggplot2)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

data_process <- function(name){
  #name = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_chr10_with_chr_pos.txt"
  count = read.delim(name, header = FALSE, sep=' ')
  count_sum <- rowSums(count[,3:1492])
  SNP <- paste0(count$V1, "_", count$V2)
  count_df <- data.frame(count$V1, count$V2, count_sum, SNP)
  outlier <- read.delim("snps_paralogs_fdr05.txt", header = FALSE, sep=' ')
  outlier$SNP <- paste0(outlier$V1, "_", outlier$V2)
  
  count_df$type[count_df$SNP %in% outlier$SNP] = "paralogs"
  count_df$type[!(count_df$SNP %in% outlier$SNP)] = "rest_SNPs"
  
  df1 <- data_summary(count_df, varname="count_sum", 
                      groupnames=c("type"))
  return(df1)
}

d1 <- data_process("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_chr9_with_chr_pos.txt")
d2 <- data_process("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_chr10_with_chr_pos.txt")

d1$chr = "chr9"
d2$chr = "chr10"

df <- rbind(d1, d2)

p<- ggplot(df, aes(x=type, y=count_sum, fill=chr)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_errorbar(aes(ymin=count_sum-sd, ymax=count_sum+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic()+
  ylab("Read count across 1490 sequenced samples")+
  xlab(NULL)+
  theme(text=element_text(family="Times New Roman", face="bold", size=14, colour="black"))

  
