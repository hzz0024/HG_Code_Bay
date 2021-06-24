library("ggplot2")
library("cowplot")
library("scales")
library("MASS")
source("individual_pca_functions.R")
file = 'pop_info.txt' 
all_label = read.delim(file, header = TRUE, sep='\t')

file1 = 'Del19_info.txt' 
Del19_label = read.delim(file1, header = TRUE, sep='\t')

file2 = 'Del20_info.txt' 
Del20_label = read.delim(file2, header = TRUE, sep='\t')

file3 = 'challenge_info.txt' 
challenge_label = read.delim(file3, header = TRUE, sep='\t')

PCA(cov_matrix = "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.covMat",
    ind_label = all_label$ind,
    pop_label = all_label$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    #index_exclude=(1:234)
)

PCoA(dist_matrix = "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)

PCA(cov_matrix = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_shared_sites.covMat",
    ind_label = all_label$ind,
    pop_label = all_label$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    #index_exclude=(1:234)
)

PCoA(dist_matrix = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_shared_sites.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)

PCoA(dist_matrix = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_shared_sites.single.read.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)

PCA(cov_matrix = "All_0p7x_minq20_minmq30_CV30_masked_noinvers_shared_sites.covMat",
    ind_label = all_label$ind,
    pop_label = all_label$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    #index_exclude=(1:234)
    )

PCoA(dist_matrix = "All_0p7x_minq20_minmq30_CV30_masked_noinvers_shared_sites.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)


PCA(cov_matrix = "All_1x_minq20_minmq30_CV30_masked_noinvers_shared_sites.covMat",
    ind_label = all_label$ind,
    pop_label = all_label$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = F,
    show.ellipse = F,
    show.line = F,
    alpha = 0.1,
    #index_exclude=(1:234)
)

PCoA(dist_matrix = "All_1x_minq20_minmq30_CV30_masked_noinvers_shared_sites.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)

### doIBS 2
PCoA(dist_matrix = "Challenge_1x_minq20_minmq30_CV30_masked_noinvers_shared_sites.ibsMat",
     ind_label = challenge_label$ind,
     pop_label = challenge_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = F,
     show.ellipse = F,
     show.line = F,
     alpha = 0.1,
     #index_exclude=(1:234)
)

