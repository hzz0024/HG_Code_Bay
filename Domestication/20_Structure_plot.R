#############################
# Plotting STRUCTURE results 
#############################

### Goal: plot Q tables from STRUCTURE analyses across 
### cross all populations (best K = 2) above separate analysis 
### subsets  within the Gulf of Mexico (K = 3) and the Atlantic Coast (K = 3) 
### populations using 10K random neutral SNPs.

### Set working directory---------------------------------------------------

setwd('./population_structure/')

### Install and load libraries/dependencies---------------------------------

#devtools::install_github('royfrancis/pophelper')
#install.packages('viridis')
library(pophelper)
library(viridis)
library(gridExtra)
library("wesanderson")
### Read in data------------------------------------------------------------
setwd("~/Dropbox/Mac/Documents/HG/Domestication/07_structure/structure_plot")

k7 <- readQ('./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_10K_structure_output_K7_run_1_f',filetype='auto')
k9 <- readQ('./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_10K_structure_output_K9_run_1_f',filetype='auto')
k14 <- readQ('./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_10K_structure_output_K14_run_1_f',filetype='auto')

# Sample metadata
inds_grps <- read.table('sampleID_pop_n509.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k7[[1]]) <- inds_grps$V1
rownames(k9[[1]]) <- inds_grps$V1
rownames(k14[[1]]) <- inds_grps$V1

### Plot CV error------------------------------------------------------------

# par(mfrow=c(1,1))
# plot(CVerror$K,CVerror$CV_error,pch=20,xlab='K',ylab='CV error')
# low <- min(CVerror$CV_error)
# abline(h=low,col='navy',lty=2,lwd=0.5)

### Plot ADMIXTURE plots-----------------------------------------------------

# set population order
pop_order <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")

# col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd",                    #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
#                    "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
                                  # blue                        # yellow
col_gradient <- c(  "#FA812F", "#849cc1",  "#1D92BD", "#8ad5d9","#fec155", "#0A2C86", "#FA4032",
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
# Groups individuals by subspecies abbreviation

p7 <- plotQ(k7, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p7$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c(  "#FA812F", "#849cc1",  "#1D92BD", "#8ad5d9","#FAAF08", "#0A2C86", "#FA4032",
                    "#f9476b", "#fec155","#fddae1", "#cf7fbc",  "#e2b2d6", "#e1bb94", "#e1bb94", "#fbd0a5", "#b58383")

p9 <- plotQ(k9, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p9$plot[[1]], nrow=1, widths=c(20,1))


col_gradient <- c(  "#FA812F", "#849cc1",  "#1D92BD", "#8ad5d9","#FAAF08", "#fec155", "#FA4032",
                    "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#0A2C86", "#e1bb94", "#fbd0a5", "#b58383")

p14 <- plotQ(k14, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.5, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(export)

jpeg("Structure_509.jpg", width = 20, height = 10, units = 'in', res = 300)
plot_grid(p7$plot[[1]],p9$plot[[1]],p14$plot[[1]], nrow=3, rel_heights = c(1/4, 1/4, 1/3))
dev.off()


#grid.arrange(p7$plot[[1]],p9$plot[[1]],p14$plot[[1]], nrow=3, rel_heights = c(1/3, 1/3, 1/3))


