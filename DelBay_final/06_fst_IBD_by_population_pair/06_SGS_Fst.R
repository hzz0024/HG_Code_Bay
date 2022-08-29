# this is a script to extract outliers for pairwise Fst calculation
#setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/06_pairwise_fst_IBD/02_SGS_outlier_fst")
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#usage: Rscript 06_SGS_Fst.R test.txt 19_SGS_HC_NB_FDR_outlier.list

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("alpha/beta file and outlier file must be supplied .n", call.=FALSE)
} 

#alpha_beta_file = "test.txt"
#outlier_list = "19_SGS_HC_NB_FDR_outlier.list"
Fst_cal <- function(alpha_beta_file, outlier_list){
  dt <- read.delim(alpha_beta_file, header = FALSE, sep='\t')
  head(dt)
  outlier_dt <- read.delim(outlier_list, header = TRUE, sep='\t')
  head(outlier_dt)
  extract_snps <- dt[which(paste0(dt$V1, "_", dt$V2) %in% outlier_dt$id), ] # for outlier fst
  #extract_snps <- dt[which(!(paste0(dt$V1, "_", dt$V2) %in% outlier_dt$id)), ] # for neutral fst
  fst_ <- sum(extract_snps$V3)/sum(extract_snps$V4)
  write.table(fst_, file = paste0(strsplit(alpha_beta_file, '[.]')[[1]][1], ".outlier.fst.txt"), sep = "\t", quote = FALSE,row.names = F, col.names = F)
}

Fst_cal(alpha_beta_file = args[1], outlier_list = args[2])


