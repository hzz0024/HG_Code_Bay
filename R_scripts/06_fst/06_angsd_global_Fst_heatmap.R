##################################################
########## global fst pairwise contrasts  ########
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/Global_pairwise_fst")

install_github("jokergoo/ComplexHeatmap")
install.packages("dendsort")
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)

data_process <- function(headname){
  # headname = "CHR19_REF19"
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.average_fst.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  #return(DT$V1) # unweighted
  return(DT$V2) # weighted, ratio between the sum of As and the sum of B
}

#data_process("HCVA_NEH")
colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}
#colnames= swap(colnames, 10, 12)
#colnames= swap(colnames, 11, 13)

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.average_fst.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname)
    else
      headname = paste0(colnames[j],'_',colnames[i])
      dat[i,j] = data_process(headname)
  }
}
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

write.table(dat2, file = "global_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0003, 0.0006), c("white", "yellow", "red")),
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

graph2ppt(file="heatmap_final",width=10,height=6)

##################################################
########## only outlier Fst contrast      ########
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/outlier_pairwise_fst")

data_process <- function(headname, outlier_list){
  #headname = "ARN_HC"
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = TRUE, sep='\t')
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
#colnames= swap(colnames, 10, 12)
#colnames= swap(colnames, 11, 13)

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "Fisher_14_outlier.txt")
    else
      headname = paste0(colnames[j],'_',colnames[i])
    dat[i,j] = data_process(headname, "Fisher_14_outlier.txt")
  }
}
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

write.table(dat2, file = "14_outlier_fst.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.05, 0.1), c("white", "yellow", "red")),
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

graph2ppt(file="outlier_heatmap_final",width=10,height=6)



##################################################
########## only pruned Fst contrast       ########
##################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/pruned_pairwise_fst/")

data_process <- function(headname, outlier_list){
  #headname = "ARN_HC"
  name = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  DT$SNP = paste0(DT$V1,'_',DT$V2)
  # read the outlier list
  OT = read.delim(outlier_list, header = FALSE, sep='\t')
  DT_outlier = DT[which(DT$SNP %in% OT$V5),]
  return(sum(DT_outlier$V3)/sum(DT_outlier$V4)) # weighted
}

colnames <- (read.delim("pop_order.txt", header = FALSE, sep='\t'))$V1
swap <- function(colnames, i, j){
  tmp = colnames[i]
  colnames[i] = colnames[j]
  colnames[j] = tmp
  return(colnames)
}
#colnames= swap(colnames, 10, 12)
#colnames= swap(colnames, 11, 13)

idx = 0
dat = matrix(0, nrow=length(colnames), ncol=length(colnames))
for(i in seq(length(colnames)-1)){
  for(j in seq(i+1, length(colnames))){
    idx = idx + 1
    print(idx)
    headname = paste0(colnames[i],'_',colnames[j])
    filename = paste0(headname, "_all_minq20_minmq30_CV30_masked_fold.alpha_beta.txt")
    if(file.exists(filename))
      dat[i,j] = data_process(headname, "Del19_pruned_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
    else
      headname = paste0(colnames[j],'_',colnames[i])
    dat[i,j] = data_process(headname, "Del19_pruned_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
  }
}
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

write.table(dat2, file = "purned_fst.csv", sep = ",", quote = FALSE,
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

graph2ppt(file="pruned_heatmap_final",width=10,height=6)


##################################################
########## only SGS candidates HC-NB      ########
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
