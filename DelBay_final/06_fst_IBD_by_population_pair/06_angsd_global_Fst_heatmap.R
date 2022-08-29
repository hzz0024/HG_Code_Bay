##################################################
########## print out the pairwise contrasts ######
##################################################

x1 <- c("HC_18", "ARN_18", "COH_18", "SR_18", "NB_18")
x2 <- c("HC_19", "ARN_19", "COH_19", "SR_19", "NB_19") 
x3 <- c("HC_21", "ARN_21", "COH_21", "SR_21", "BS_21", "BEN_21", "NAN_21", "NB_21")
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
    my_vec <- c(my_vec, paste("pop", comb[i,][1], comb[i,][2], sep="_"))    # Appending new value to vector
    #my_vec <- c(my_vec, paste(comb[i,][1], comb[i,][2], sep="_"))    # Appending new value to vector
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
#setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/06_Fst/Global_pairwise_fst")
rm(list=ls())
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/08_prune_IBD/WILD21/")

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
  #name = paste0(headname, "_sfs_all_sites_fold.average_fst.txt")
  DT = read.delim(name, header = FALSE, sep='\t')
  #return(DT$V1) # unweighted
  return(DT$V2) # weighted, ratio between the sum of As and the sum of B
}

#data_process("HCVA_NEH")
colnames <- (read.delim("pop_order_WILD21.txt", header = FALSE, sep='\t'))$V1

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

write.table(dat2, file = "WILD18.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('ggplot2')

Heatmap(dat2, name = "Global Fst", row_dend_reorder = FALSE, column_dend_reorder=FALSE, column_title = "Population", 
        row_order = order(as.numeric(gsub("row", "", rownames(dat2)))),
        column_order = order(as.numeric(gsub("column", "", colnames(dat2)))),
        col = colorRamp2(c(0, 0.0005, 0.001), c("white", "yellow", "red")),
        #col = colorRamp2(c(0, 0.0005, 0.001), c("white", "yellow", "red")),
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

graph2ppt(file="WILD21",width=10,height=6)

