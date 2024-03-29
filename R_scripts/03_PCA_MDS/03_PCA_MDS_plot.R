library("ggplot2")
library("cowplot")
library("scales")
library("MASS")
library(export)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/03_global_PCA_MDS/")
source("individual_pca_functions.R")
# for color palettes see https://r-charts.com/color-palettes/#discrete > paletteer_d("ggthemes::Hue_Circle")
file = 'Del20_101_info.txt' 
all_label = read.delim(file, header = TRUE, sep='\t')

PCoA(dist_matrix = "Del20_challenge_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers_nodownsampling.ibsMat",
     ind_label = all_label$ind,
     pop_label = all_label$pop,
     x_axis = 1,
     y_axis = 2,
     k = 2,
     show.point = T,
     show.label = T,
     show.ellipse = T,
     show.line = F,
     alpha = 0.05,
)

PCA(cov_matrix = "All_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers_shared_sites.covMat",
    ind_label = all_label$ind,
    pop_label = all_label$pop,
    x_axis = 1,
    y_axis = 2,
    show.point = T,
    show.label = T,
    show.ellipse = T,
    show.line = F,
    alpha = 0.1,
    #index_exclude=(1:234)
)

file = 'CLP_15.txt' 
CLP = read.delim(file, header = TRUE, sep='\t')

PCoA(dist_matrix = "CLP_15_minmapq30_minq20_CV30_masked_noinvers_shared_sites.ibsMat",
     ind_label = CLP$ind,
     pop_label = CLP$pop,
     x_axis = 1,
     y_axis = 2,
     k = 3,
     show.point = T,
     show.label = T,
     show.ellipse = T,
     show.line = T,
     alpha = 0.05,
     #index_exclude=(1:38)
)
graph2ppt(file="All_MDS.pptx", width=4, height=3)

file = 'All_92_info.txt' 
all_label = read.delim(file, header = TRUE, sep='\t')

graph2ppt(file="All_PCA.pptx", width=4, height=3)
############################################
##### alternative way of MDS plotting ######
############################################
library(ggplot2)
library(ggrepel)
#challenge
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ff0000", "#ff7878", "#D55E00", "#CC79A7", "#8FDDE7")

name <- "All_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_noinvers_shared_sites.ibsMat"
m <- as.matrix(read.table(name))
mds <- cmdscale(as.dist(m))
anno <- read.table('All_432_list.txt', sep="\t", header=T)
mds <- as.data.frame(mds)
mds$Population <- anno$CLUSTER

# define the range index for labelling
#idx_in_range0 = mds$V1>-0.0 & mds$V1<0 & mds$V2 >-0.0 & mds$V2 < -0.0
#idx_in_range1 = mds$V1>0.0 
#idx_in_range = idx_in_range0 | idx_in_range1
idx_in_range = mds$V1> 0.3
print(idx_in_range)

#plot(mds,lwd=2,ylab="Dimension 2",xlab="Dimension 1",main="Multidimensional Scaling Plot",col=rep(1:5,each=1))
ggplot() + geom_point(data=mds, aes_string(x='V1', y='V2', color="Population"))+ 
    ggtitle('Multidimensional Scaling')+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Dimension 1") + ylab("Dimension 2") +
    scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB", "CHR19","REF19","CHR20","REF20"))+
    geom_text_repel(data=mds[idx_in_range,], aes_string(x='V1', y='V2', color="Population"), label=anno$IID[idx_in_range], size=2, nudge_x = 0.001, nudge_y = 0.001)
ggsave("All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites.jpg")

# # wild
# cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ff0000", "#ff7878", "#D55E00", "#CC79A7")
# 
# name <- "wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.ibsMat"
# m <- as.matrix(read.table(name))
# mds <- cmdscale(as.dist(m))
# anno <- read.table('wild_235.list', sep="\t", header=T)
# mds <- as.data.frame(mds)
# mds$population <- anno$CLUSTER
# 
# # define the range index for labelling
# idx_in_range0 = mds$V1>0.02 & mds$V1<0.06 & mds$V2 >-0.04 & mds$V2 < 0.03
# idx_in_range1 = mds$V2>0.03 
# idx_in_range = idx_in_range0 | idx_in_range1
# print(idx_in_range)
# 
# #plot(mds,lwd=2,ylab="Dimension 2",xlab="Dimension 1",main="Multidimensional Scaling Plot",col=rep(1:5,each=1))
# ggplot() + geom_point(data=mds, aes_string(x='V1', y='V2', color="population"))+ 
#   ggtitle('Multidimensional Scaling')+
#   theme(plot.title = element_text(hjust = 0.5))+
#   xlab("Dimension 1") + ylab("Dimension 2") +
#   scale_colour_manual(values=cbPalette, breaks=c("HC","ARN","COH","SR","NB"))+
#   geom_text_repel(data=mds[idx_in_range,], aes_string(x='V1', y='V2', color="population"), label=anno$IID[idx_in_range], size=2, nudge_x = 0.001, nudge_y = 0.001)

ggsave("wild_MDS.jpg")
