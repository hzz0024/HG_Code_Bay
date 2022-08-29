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

devtools::install_github('royfrancis/pophelper')
install.packages('viridis')
library(pophelper)
library(viridis)
library(gridExtra)
### Read in data------------------------------------------------------------
setwd("~/Dropbox/Mac/Documents/HG/CVreseq/STRUCTURE/Sturcutre_plot")
# CV error
CVerror <- read.table('./CV_error.txt',header=T)

# Q tables
k2 <- readQ('./all_n_78_random10k_structure_output_K2_run_1_f',filetype='auto')
k3 <- readQ('./cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.3.Q',filetype='auto')
k4 <- readQ('./cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.4.Q',filetype='auto')

# Sample metadata
inds_grps <- read.table('sampleID_pop_n78.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k2[[1]]) <- inds_grps$V1
rownames(k3[[1]]) <- inds_grps$V1
rownames(k4[[1]]) <- inds_grps$V1
rownames(k5[[1]]) <- inds_grps$V1

### Plot CV error------------------------------------------------------------

# par(mfrow=c(1,1))
# plot(CVerror$K,CVerror$CV_error,pch=20,xlab='K',ylab='CV error')
# low <- min(CVerror$CV_error)
# abline(h=low,col='navy',lty=2,lwd=0.5)

### Plot ADMIXTURE plots-----------------------------------------------------

# set population order
pop_order <- c("CL","SL","LOLA","OBOYS","HC","CS","CLP","HC_VA","HI","SM","DEBY","UMFS","NEH")

# Groups individuals by subspecies abbreviation

p1 <- plotQ(k2, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=c("#D09683","#92AAC7","lightblue","purple","orange"),
            barsize = 0.9, barbordercolour="black",barbordersize=0.6,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p1$plot[[1]], nrow=1, widths=c(10,1))

k3 <- readQ('./GM_n_24_random10k_structure_output_K3_run_1_f',filetype='auto')
pop_order <- c("CL","SL","LOLA","OBOYS")

# Sample metadata
inds_grps <- read.table('sampleID_pop_n24.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k3[[1]]) <- inds_grps$V1
p2 <- plotQ(k3, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=c("#EE693F","#D09683","#752A07"),
            barsize = 1, barbordercolour="black",barbordersize=0.6,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p2$plot[[1]], nrow=1, widths=c(10,1))


k3 <- readQ('./AC_n_48_random10k_structure_output_K3_run_1_f',filetype='auto')
pop_order <- c("HC","CS","CLP","HC_VA","HI","SM","UMFS","NEH")

# Sample metadata
inds_grps <- read.table('sampleID_pop_n48.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k3[[1]]) <- inds_grps$V1
p2 <- plotQ(k3, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=c("#4D648D","#92AAC7","#1E1F26"),
            barsize = 1, barbordercolour="black",barbordersize=0.6,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p2$plot[[1]], nrow=1, widths=c(10,1))
