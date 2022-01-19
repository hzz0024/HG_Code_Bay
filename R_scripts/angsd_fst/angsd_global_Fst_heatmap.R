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
  return(DT$V1) # weighted
  #return(DT$V2) # unweighted 
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
