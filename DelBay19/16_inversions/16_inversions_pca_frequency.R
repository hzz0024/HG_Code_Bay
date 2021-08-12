library("ggplot2")
library("cowplot")
library("scales")
library("MASS")
setwd("~/Documents/HG/DelBay19_adult/16_inversions")
source("individual_pca_functions.R")

file = 'HC_info.txt' 
HC = read.delim(file, header = TRUE, sep='\t')

file = 'ARN_info.txt' 
ARN = read.delim(file, header = TRUE, sep='\t')

file = 'COH_info.txt' 
COH = read.delim(file, header = TRUE, sep='\t')

file = 'SR_info.txt' 
SR = read.delim(file, header = TRUE, sep='\t')

file = 'NB_info.txt' 
NB = read.delim(file, header = TRUE, sep='\t')

file = 'REF19_info.txt' 
REF19 = read.delim(file, header = TRUE, sep='\t')

file = 'CHR19_info.txt' 
CHR19 = read.delim(file, header = TRUE, sep='\t')


########### HC ###########
source("individual_pca_functions.R")
jpeg("chr1_HC.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_HC_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = HC$ind,
    pop_label = HC$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="HC chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_HC.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_HC_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = HC$ind,
    pop_label = HC$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="HC chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_HC.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_HC_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = HC$ind,
    pop_label = HC$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="HC chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### ARN ###########
source("individual_pca_functions.R")
jpeg("chr1_ARN.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_ARN_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = ARN$ind,
    pop_label = ARN$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="ARN chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_ARN.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_ARN_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = ARN$ind,
    pop_label = ARN$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="ARN chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_ARN.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_ARN_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = ARN$ind,
    pop_label = ARN$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="ARN chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### COH ###########
source("individual_pca_functions.R")
jpeg("chr1_COH.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_COH_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = COH$ind,
    pop_label = COH$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="COH chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_COH.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_COH_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = COH$ind,
    pop_label = COH$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="COH chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_COH.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_COH_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = COH$ind,
    pop_label = COH$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="COH chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### SR ###########
source("individual_pca_functions.R")
jpeg("chr1_SR.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_SR_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = SR$ind,
    pop_label = SR$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="SR chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_SR.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_SR_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = SR$ind,
    pop_label = SR$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="SR chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_SR.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_SR_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = SR$ind,
    pop_label = SR$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="SR chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### NB ###########
source("individual_pca_functions.R")
jpeg("chr1_NB.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_NB_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = NB$ind,
    pop_label = NB$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="NB chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_NB.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_NB_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = NB$ind,
    pop_label = NB$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="NB chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_NB.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_NB_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = NB$ind,
    pop_label = NB$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="NB chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### CHR19 ###########
source("individual_pca_functions.R")
jpeg("chr1_CHR19.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_CHR19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = CHR19$ind,
    pop_label = CHR19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="CHR19 chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_CHR19.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_CHR19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = CHR19$ind,
    pop_label = CHR19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="CHR19 chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_CHR19.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_CHR19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = CHR19$ind,
    pop_label = CHR19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="CHR19 chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### REF19 ###########
source("individual_pca_functions.R")
jpeg("chr1_REF19.jpg", width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1,1))
#par(mar=c(4,8,1,6))
PCA(cov_matrix = "chr1_REF19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = REF19$ind,
    pop_label = REF19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="REF19 chromosome 1 inversions (40630000-42400000)"
)
dev.off()
jpeg("chr5_REF19.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr5_REF19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = REF19$ind,
    pop_label = REF19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="REF19 chromosome 5 inversions (62190000-79000000)"
)
dev.off()
jpeg("chr6_REF19.jpg", width = 6, height = 4, units = 'in', res = 300)
PCA(cov_matrix = "chr6_REF19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.covMat",
    ind_label = REF19$ind,
    pop_label = REF19$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    main="REF19 chromosome 6 inversions (30500000-43900000)"
)
dev.off()

########### calculate the R2 between frequency and salinity index ###########




