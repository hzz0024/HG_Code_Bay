rm(list=ls())

library(dplyr)
library(VennDiagram)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/")


########################
# collect RDA outliers #
########################

rda_outlier <- function(predictor){
  #load outliers tables
  #predictor = "SA"
  outlier_temp_rda<-read.table("03_rda/2.5sd/info_pop_env_rda.txt_candidate_SNP_2.5_sd.txt", header=T)
  head(outlier_temp_rda)
  nRDA<-dim(outlier_temp_rda)[1]
  nRDA #how many outliers?
  target_outlier <- outlier_temp_rda[which(outlier_temp_rda$predictor == predictor),]
  print(paste0("Number of rda outliers for ", predictor, " is ", dim(target_outlier)[1])) # count how many SNPs with predictor factor
  return(target_outlier$snp)
  #assign(paste0(predictor,".rda.outliers"), target_outlier, envir = globalenv())
}

rda_outlier("SA")
rda_outlier("CD5")
rda_outlier("CD7")
rda_outlier("CD9")
rda_outlier("CD11")
rda_outlier("MAX10")

#########################
# collect lfmm outliers #
#########################

lfmm_outlier <- function(fname) {
  #fname = "CD11_k1" 
  outlier_temp_lfmm <- read.delim(paste0("04_lfmm/all_env/", fname, ".FDR.outlier.txt"), header = F, sep="\t")
  colnames(outlier_temp_lfmm) <- c("chr", "pos", "fdr", "id")
  outlier_temp_lfmm$pos <- as.numeric(outlier_temp_lfmm$pos * 1e7)
  # count how many SNPs with FDR < fdr_cut
  print(paste0("Number of lfmm outliers for ", fname, " is ", length(outlier_temp_lfmm$id))) # count how many SNPs with predictor factor
  # format for common snps
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_temp_lfmm$chr[outlier_temp_lfmm$chr==i] = chr_str_list[i]
  outlier_temp_lfmm$snp_id <- paste0(outlier_temp_lfmm$chr, "_", outlier_temp_lfmm$pos)
  return(outlier_temp_lfmm$snp_id)
  #assign(paste0(fname,".lfmm.outliers"), outlier_temp_lfmm, envir = globalenv())
}

lfmm_outlier("SA_k1")
lfmm_outlier("CD5_k1")
lfmm_outlier("CD7_k1")
lfmm_outlier("CD9_k1")
lfmm_outlier("CD11_k1")
lfmm_outlier("MAX10_k1")

############################
# collect Baypass outliers #
############################

baypass_outlier <- function(predictor){
  #load outliers tables
  #predictor = "SA"
  outlier_temp_bp<-read.table(paste0("02_baypass/",predictor, ".baypass.FDR.outlier.txt"), header=F)
  head(outlier_temp_bp)
  
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_temp_bp$V1[outlier_temp_bp$V1==i] = chr_str_list[i]
  outlier_temp_bp$snp_id <- paste0(outlier_temp_bp$V1, "_", outlier_temp_bp$V2*1e+7)
  nBP<-dim(outlier_temp_bp)[1] #how many outliers?
  print(paste0("Number of rda outliers for ", predictor, " is ", nBP)) # count how many SNPs with predictor factor
  return(outlier_temp_bp$snp_id)
}

baypass_outlier("SA")
baypass_outlier("CD5")
baypass_outlier("CD7")
baypass_outlier("CD9")
baypass_outlier("CD11")
baypass_outlier("MAX10")
######################
# output shared snps #
######################

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

for(pop in c("SA","CD5","CD7","CD9","CD11","MAX10")){
  shared <- intersect(rda_outlier(pop), lfmm_outlier(paste0(pop, "_k1")))
  chr1 <- function(shared){
    unlist(strsplit(shared, "_"))[1]
  }
  chr2 <- function(shared){
    unlist(strsplit(shared, "_"))[2]
  }
  pos <- function(shared){
    unlist(strsplit(shared, "_"))[3]
  }
  shared_df <- cbind(paste0(sapply(shared, chr1),"_", sapply(shared, chr2)), sapply(shared, pos), sapply(shared, pos))
  shared_df[,2] <- as.numeric(shared_df[,2]) 
  shared_df[,3] <- as.numeric(shared_df[,3]) 
  #shared_df[,2] <- as.numeric(shared_df[,2]) - 5000 # for 10K range of the snps
  #shared_df[,3] <- as.numeric(shared_df[,3]) + 5000 # for 10K range of the snps
  
  write.table(shared_df, paste0("./06_shared/", pop,"_shared_lfmm_rda.bed"), row.names = FALSE,  col.names = FALSE, sep="\t", quote = FALSE)
  x <- list(A = rda_outlier(pop), 
            B = lfmm_outlier(paste0(pop, "_k1")))
  display_venn(x,
               lwd = 1,
               lty = 'blank',
               category.names = c( paste0(pop, "_rda_outlier") , paste0(pop, "_lfmm_outlier")),
               fill = c("#999999", "#E69F00"),
               cat.cex = 1,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.dist = c(0.12, 0.27))
  print(paste0("Number of shared outlier between rda and lfmm is ", length(shared)))
}



#############################################
# output shared snps among three approaches #
#############################################

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

for(pop in c("SA","CD5","CD7","CD9","CD11","MAX10")){
  shared <- intersect(intersect(rda_outlier(pop), lfmm_outlier(paste0(pop, "_k1"))), baypass_outlier(pop))
  chr1 <- function(shared){
    unlist(strsplit(shared, "_"))[1]
  }
  chr2 <- function(shared){
    unlist(strsplit(shared, "_"))[2]
  }
  pos <- function(shared){
    unlist(strsplit(shared, "_"))[3]
  }
  if(length(shared) > 0){
    shared_df <- cbind(paste0(sapply(shared, chr1),"_", sapply(shared, chr2)), sapply(shared, pos), sapply(shared, pos))
    shared_df[,2] <- as.numeric(shared_df[,2]) 
    shared_df[,3] <- as.numeric(shared_df[,3]) 
    #shared_df[,2] <- as.numeric(shared_df[,2]) - 5000 # for 10K range of the snps
    #shared_df[,3] <- as.numeric(shared_df[,3]) + 5000 # for 10K range of the snps
    
    write.table(shared_df, paste0("./06_shared/", pop,"_shared_lfmm_rda_baypass.bed"), row.names = FALSE,  col.names = FALSE, sep="\t", quote = FALSE)
    x <- list(A = rda_outlier(pop), 
              B = lfmm_outlier(paste0(pop, "_k1")),
              C = baypass_outlier(pop))
    display_venn(x,
                 lwd = 1,
                 lty = 'blank',
                 category.names = c( paste0(pop, "_rda_outlier") , paste0(pop, "_lfmm_outlier"), paste0(pop, "_baypass_outlier")),
                 fill = c("#999999", "#E69F00", "red"),
                 cat.cex = 1,
                 cat.fontface = "bold",
                 cat.default.pos = "outer",
                 cat.dist = c(0.12, 0.12, 0.12))
    print(paste0("Number of shared outliers among rda, lfmm, and baypass is ", length(shared)))
  }
}

##############################
# random test on shared SNPs #
##############################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/05_random_test")
snp_list <- read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", header = F, sep="\t")
snps <- paste0(snp_list$V1, "_", snp_list$V2)

rep=10000
cnt1 = 1694
cnt2 = 523
target =45
cnt <- c()
for (i in 1:rep){
  d1 <- sample(snps, cnt1, replace=FALSE)
  d2 <- sample(snps, cnt2, replace=FALSE)
  share_cnt <- length(intersect(d1,d2))
  cnt <- c(cnt, share_cnt)
}
hist(cnt)
sum(cnt > target)/rep
  






