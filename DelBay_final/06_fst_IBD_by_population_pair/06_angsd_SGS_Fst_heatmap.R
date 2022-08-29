##################################################
########## print out the pairwise contrasts ######
##################################################

x1 <- c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18")
x2 <- c("HC_19", "ARN_19", "COH_19", "SR_19", "NB_19") 
x3 <- c("HC_21", "ARN_21", "BEN_21", "BS_21", "COH_21", "SR_21", "NAN_21", "NB_21")
x4 <- c("HC_21spat", "ARN_21spat", "NAN_21spat")
x5 <- c("HC_18", "HC_19", "HC_21", "HC_21spat")
x6 <- c("ARN_18", "ARN_19", "ARN_21", "ARN_21spat")
x7 <- c("NAN_21", "NAN_21spat")
x8 <- c("Sur_19", "Ref_19")
x9 <- c("Sur_20", "Ref_20")
x10 <- c("A_Sur20", "B_Sur20", "Ref_20", "Ref_19")

pop_vector <- function (x) {
  n <- 2
  comb <- t(combn(x,n))
  my_vec <- c()
  for(i in 1:dim(comb)[1]) {
    #my_vec <- c(my_vec, paste(comb[i,][1], comb[i,][2], sep="_"))    # Appending new value to vector
    my_vec <- c(my_vec, paste("pop", comb[i,][1], comb[i,][2], sep="_"))    # Appending new value to vector
  }
  print(noquote(my_vec))
} 

pop_vector(x1)
pop_vector(x2)
pop_vector(x3)
pop_vector(x4)
pop_vector(x5)
pop_vector(x6)
pop_vector(x7)
pop_vector(x8)
pop_vector(x9)
pop_vector(x10)

##################################################
########## global fst pairwise contrasts  ########
##################################################
rm(list=ls())
#install_github("jokergoo/ComplexHeatmap")
#install.packages("dendsort")
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)

##################################################
########## only SGS candidates HC-NB      ########
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/08_fst_IBD/SGS_outlier_Fst")

data_process <- function(headname, outlier_list){
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = TRUE, sep='\t')
  message("Number of target list has ", dim(OT)[1], " SNPs")
  DT_outlier = DT[which(DT$SNP %in% OT$id),]
  return(sum(DT_outlier$V3)/sum(DT_outlier$V4)) # weighted
}

colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "HC_NB_FDR_outlier.list")
    else
      headname = paste0(colnames[j],'_',colnames[i])
      dat[i,j] = data_process(headname, "HC_NB_FDR_outlier.list")
  }
}
dim(dat)
colnames(dat) = colnames
rownames(dat) = colnames
dat1 = dat
dat1[dat1==0]=NA

dat2 = dat
for(i in seq(length(colnames)))
  dat2[i, i] = 0 
for(i in seq(length(colnames)-1))
  for(j in seq(i+1, length(colnames)))
    dat2[j,i] = dat2[i,j]

write.table(dat2, file = "SGS_HC_NB_outlier_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0002, 0.0005), c("white", "yellow", "red")),
        #row_dend_height = unit(2, "cm"),
        #km = 2,
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", dat2[i, j]), x, y, gp = gpar(fontsize = 15))
        })

graph2ppt(file="SGS_candidate_HC_NB_heatmap_final",width=10,height=6)

##################################################
########## only SGS candidates CHR19-REF19 #######
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/outlier_pairwise_fst/")

data_process <- function(headname, outlier_list){
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = TRUE, sep='\t')
  message("Number of target list has ", dim(OT)[1], " SNPs")
  DT_outlier = DT[which(DT$SNP %in% OT$id),]
  return(sum(DT_outlier$V3)/sum(DT_outlier$V4)) # weighted
}

colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "challenge_FDR_outlier.list")
    else
      headname = paste0(colnames[j],'_',colnames[i])
    dat[i,j] = data_process(headname, "challenge_FDR_outlier.list")
  }
}
dim(dat)
colnames(dat) = colnames
rownames(dat) = colnames
dat1 = dat
dat1[dat1==0]=NA

dat2 = dat
for(i in seq(length(colnames)))
  dat2[i, i] = 0 
for(i in seq(length(colnames)-1))
  for(j in seq(i+1, length(colnames)))
    dat2[j,i] = dat2[i,j]

write.table(dat2, file = "SGS_CHR19_REF19_outlier_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0002, 0.0005), c("white", "yellow", "red")),
        #row_dend_height = unit(2, "cm"),
        #km = 2,
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", dat2[i, j]), x, y, gp = gpar(fontsize = 15))
        })

graph2ppt(file="SGS_candidate_CHR19_REF19_heatmap_final",width=10,height=6)


##########################################################
########## both SGS candidates HC-NB + CHR19-REF19 #######
##########################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/outlier_pairwise_fst/")

data_process <- function(headname, outlier_list){
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = TRUE, sep='\t')
  message("Number of target list has ", dim(OT)[1], " SNPs")
  DT_outlier = DT[which(DT$SNP %in% OT$id),]
  return(sum(DT_outlier$V3)/sum(DT_outlier$V4)) # weighted
}

colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "both_SGS_FDR_outlier.list")
    else
      headname = paste0(colnames[j],'_',colnames[i])
    dat[i,j] = data_process(headname, "both_SGS_FDR_outlier.list")
  }
}
dim(dat)
colnames(dat) = colnames
rownames(dat) = colnames
dat1 = dat
dat1[dat1==0]=NA

dat2 = dat
for(i in seq(length(colnames)))
  dat2[i, i] = 0 
for(i in seq(length(colnames)-1))
  for(j in seq(i+1, length(colnames)))
    dat2[j,i] = dat2[i,j]

write.table(dat2, file = "SGS_both_outlier_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0002, 0.0005), c("white", "yellow", "red")),
        #row_dend_height = unit(2, "cm"),
        #km = 2,
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", dat2[i, j]), x, y, gp = gpar(fontsize = 15))
        })

graph2ppt(file="SGS_candidate_CHR19_REF19_heatmap_final",width=10,height=6)
##################################################
########## no SGS candidates HC-NB        ########
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/no_SGS_candidates/")

data_process <- function(headname, outlier_list){
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = TRUE, sep='\t')
  DT_outlier = DT[which(!DT$SNP %in% OT$id),]
  message("Number of snp is ", dim(DT_outlier)[1], " SNPs")
  return(sum(DT_outlier$V3)/sum(DT_outlier$V4)) # weighted
}

colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "HC_NB_FDR_outlier.list")
    else
      headname = paste0(colnames[j],'_',colnames[i])
      dat[i,j] = data_process(headname, "HC_NB_FDR_outlier.list")
  }
}
dim(dat)
colnames(dat) = colnames
rownames(dat) = colnames
dat1 = dat
dat1[dat1==0]=NA

dat2 = dat
for(i in seq(length(colnames)))
  dat2[i, i] = 0 
for(i in seq(length(colnames)-1))
  for(j in seq(i+1, length(colnames)))
    dat2[j,i] = dat2[i,j]

write.table(dat2, file = "genome_no_SGS_outlier_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0002, 0.0005), c("white", "yellow", "red")),
        #row_dend_height = unit(2, "cm"),
        #km = 2,
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", dat2[i, j]), x, y, gp = gpar(fontsize = 15))
        })

graph2ppt(file="SGS_candidate_HC_NB_heatmap_final",width=10,height=6)
