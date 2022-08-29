setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Outliers/share_evalucation")
install.packages("ggVennDiagram")
library("ggVennDiagram")
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)
library(stats)

format_outlier <- function(outlier_name){
  #outlier_name="18_HC_NB_shared_outlier.txt"
  outlier_df = read.delim(outlier_name, header = TRUE, sep='\t') 
  outlier_SNP <- outlier_df$id
  return(outlier_SNP)
}

HC_NB_18 <- format_outlier("18_HC_NB_shared_outlier.txt")
HC_NB_19 <- format_outlier("19_HC_NB_shared_outlier.txt")
HC_NB_21 <- format_outlier("21_HC_NB_shared_outlier.txt")
Sur_Ref_19 <- format_outlier("19_Sur_Ref_shared_outlier.txt")
Sur_Ref_20 <- format_outlier("20_Sur_Ref_shared_outlier.txt")
CD5_df <- read.delim("CD5_shared_lfmm_rda.bed", header = FALSE, sep='\t')
GEA_CD5 <- paste0(CD5_df$V1, "_", CD5_df$V2)

x = list( HC_NB_18 = HC_NB_18,
          HC_NB_19 = HC_NB_19,
          HC_NB_21 = HC_NB_21,
          Sur_Ref_19= Sur_Ref_19,
          Sur_Ref_20 = Sur_Ref_20,
          GEA_CD5=GEA_CD5)

ggVennDiagram(x)+ scale_color_brewer(palette = "Paired")

###########################
### count for intervals ###
###########################
format_bed <- function(outlier_name, window){
  #outlier_name="18_HC_NB_shared_outlier.txt"
  outlier_df = read.delim(outlier_name, header = TRUE, sep='\t') 
  outlier_bed <- data.frame(outlier_df$chr, outlier_df$pos-window, outlier_df$pos+window)
  names(outlier_bed) <- NULL 
  #return(outlier_bed)
  write.table(outlier_bed, file = paste0(outlier_name, ".10k.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

HC_NB_18 <- format_bed("18_HC_NB_shared_outlier.txt", 500)
HC_NB_19 <- format_bed("19_HC_NB_shared_outlier.txt", 500)
HC_NB_21 <- format_bed("21_HC_NB_shared_outlier.txt", 500)
Sur_Ref_19 <- format_bed("19_Sur_Ref_shared_outlier.txt", 500)
Sur_Ref_20 <- format_bed("20_Sur_Ref_shared_outlier.txt", 500)
CD5_df <- read.delim("CD5_shared_lfmm_rda.bed", header = FALSE, sep='\t')
GEA_CD5 <- paste0(CD5_df$V1, "_", CD5_df$V2)
CD5_df_10K <- data.frame(CD5_df$V1, CD5_df$V2-window, CD5_df$V3+window)
write.table(CD5_df_10K, file = "CD5_shared_lfmm_rda.bed.10k.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bedtools = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/bedtools2/bin/bedtools"
system("bedtools intersect -a 18_HC_NB_shared_outlier.txt.10K.bed -b 19_HC_NB_shared_outlier.txt.10K.bed -wo > 18_HC_NB_shared_outlier.txt.bed_19_HC_NB_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 18_HC_NB_shared_outlier.txt.10K.bed -b 21_HC_NB_shared_outlier.txt.10K.bed -wo > 18_HC_NB_shared_outlier.txt.bed_21_HC_NB_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 18_HC_NB_shared_outlier.txt.10K.bed -b 19_Sur_Ref_shared_outlier.txt.10K.bed -wo > 18_HC_NB_shared_outlier.txt.bed_19_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 18_HC_NB_shared_outlier.txt.10K.bed -b 20_Sur_Ref_shared_outlier.txt.10K.bed -wo > 18_HC_NB_shared_outlier.txt.bed_20_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 18_HC_NB_shared_outlier.txt.10K.bed -b CD5_shared_lfmm_rda.bed.10k.bed -wo > 18_HC_NB_shared_outlier.txt.bed_CD5_shared_lfmm_rda.bed.10K.intersect", intern = TRUE)

system("bedtools intersect -a 19_HC_NB_shared_outlier.txt.10K.bed -b 21_HC_NB_shared_outlier.txt.10K.bed -wo > 19_HC_NB_shared_outlier.txt.bed_21_HC_NB_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 19_HC_NB_shared_outlier.txt.10K.bed -b 19_Sur_Ref_shared_outlier.txt.10K.bed -wo > 19_HC_NB_shared_outlier.txt.bed_19_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 19_HC_NB_shared_outlier.txt.10K.bed -b 20_Sur_Ref_shared_outlier.txt.10K.bed -wo > 19_HC_NB_shared_outlier.txt.bed_20_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 19_HC_NB_shared_outlier.txt.10K.bed -b CD5_shared_lfmm_rda.bed.10k.bed -wo > 19_HC_NB_shared_outlier.txt.bed_CD5_shared_lfmm_rda.bed.10K.intersect", intern = TRUE)

system("bedtools intersect -a 21_HC_NB_shared_outlier.txt.10K.bed -b 19_Sur_Ref_shared_outlier.txt.10K.bed -wo > 21_HC_NB_shared_outlier.txt.bed_19_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 21_HC_NB_shared_outlier.txt.10K.bed -b 20_Sur_Ref_shared_outlier.txt.10K.bed -wo > 21_HC_NB_shared_outlier.txt.bed_20_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 21_HC_NB_shared_outlier.txt.10K.bed -b CD5_shared_lfmm_rda.bed.10k.bed -wo > 21_HC_NB_shared_outlier.txt.bed_CD5_shared_lfmm_rda.bed.10K.intersect", intern = TRUE)

system("bedtools intersect -a 19_Sur_Ref_shared_outlier.txt.10K.bed -b 20_Sur_Ref_shared_outlier.txt.10K.bed -wo > 19_Sur_Ref_shared_outlier.txt.bed_20_Sur_Ref_shared_outlier.txt.bed.10K.intersect", intern = TRUE)
system("bedtools intersect -a 19_Sur_Ref_shared_outlier.txt.10K.bed -b CD5_shared_lfmm_rda.bed.10k.bed -wo > 19_Sur_Ref_shared_outlier.txt.bed_CD5_shared_lfmm_rda.bed.10K.intersect", intern = TRUE)

system("bedtools intersect -a 20_Sur_Ref_shared_outlier.txt.10K.bed -b CD5_shared_lfmm_rda.bed.10k.bed -wo > 20_Sur_Ref_shared_outlier.txt.bed_CD5_shared_lfmm_rda.bed.10K.intersect", intern = TRUE)

summary_df = read.delim("global_fst.csv", header = TRUE, sep=',') 
summary_df <-as.matrix(summary_df)

Heatmap(summary_df, name = "Shared outliers", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Contrast", 
        col = colorRamp2(c(0, 110), c("white", "red")),
        cluster_rows = FALSE,
        cluster_columns =FALSE,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", summary_df[i, j]), x, y, gp = gpar(fontsize = 8))
        })

library(dendsort)
row_dend = dendsort(hclust(dist(summary_df)))
col_dend = dendsort(hclust(dist(t(summary_df))))
Heatmap(summary_df, name = "Global Fst", cluster_rows = row_dend, cluster_columns = col_dend, column_title = "Population/Line", 
        #col = colorRamp2(c(0, 0.1, 0.2), c("white", "blue", "red")),
        col = colorRamp2(c(0, 0.2), c("white", "red")),
        #row_dend_height = unit(2, "cm"),
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", summary_df[i, j]), x, y, gp = gpar(fontsize = 8))
        })


