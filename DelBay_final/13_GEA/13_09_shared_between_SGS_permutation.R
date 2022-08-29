library(wesanderson)
rm(list=ls())

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/07_trend_plot/")
snp_df <- read.delim("by_pop_0.05_pctind0.7_maxdepth3.mafs", header = T, sep=" ")
snp_df$snps <- paste0(snp_df$chromo, "_", snp_df$position)

outlier_df <- read.delim("CD5_shared_lfmm_rda.bed", header = F, sep="\t")
outliers <- paste0(outlier_df$V1, "_", outlier_df$V2)

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/13_env_gen_association/09_shared_outliers_with_SGS&permutation")

HC_NB_18 <- read.delim("18_HC_NB_shared_outlier.txt", header = T, sep="\t")
intersect(outliers, HC_NB_18$id)

HC_NB_19 <- read.delim("19_HC_NB_shared_outlier.txt", header = T, sep="\t")
intersect(outliers, HC_NB_19$id)

HC_NB_21 <- read.delim("21_HC_NB_shared_outlier.txt", header = T, sep="\t")
intersect(outliers, HC_NB_21$id)

Sur_Ref_19 <- read.delim("19_Sur_Ref_shared_outlier.txt", header = T, sep="\t")
intersect(outliers, Sur_Ref_19$id)

Sur_Ref_21 <- read.delim("21_HC_NB_shared_outlier.txt", header = T, sep="\t")
intersect(outliers, Sur_Ref_21$id)

