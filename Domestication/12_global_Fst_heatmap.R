install_github("jokergoo/ComplexHeatmap")
install.packages("dendsort")
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/06_global_fst")
dat1 <- read.delim("global_fst.csv", header = TRUE, sep=',')
dat2 <- as.matrix(dat1)
rownames(dat2) <- colnames(dat2)
#Heatmap(dat1, name = "Fst", show_row_names = TRUE, row_order = 12:1, column_order = 12:1)
library('gplots')

Heatmap(dat2, name = "Global Fst", column_title = "Population/Line", 
        row_order = seq(1:17),
        column_order = seq(1:17),
        #col = colorRamp2(c(0, 0.1, 0.2), c("white", "blue", "red")),
        col = colorRamp2(c(0, 0.2), c("white", "red")),
        #row_dend_height = unit(2, "cm"),
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", dat2[i, j]), x, y, gp = gpar(fontsize = 8))
        })

library(dendsort)
row_dend = dendsort(hclust(dist(dat2)))
col_dend = dendsort(hclust(dist(t(dat2))))

Heatmap(dat2, name = "Global Fst", cluster_rows = row_dend, cluster_columns = col_dend, column_title = "Population/Line",
          col = colorRamp2(c(0, 0.2), c("white", "red")),
        #row_dend_height = unit(2, "cm"),
        #split = c(rep("Gulf", 4), rep("Atlantic", 9)),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        row_dend_side = "left",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", dat2[i, j]), x, y, gp = gpar(fontsize = 8))
        })


graph2ppt(file="heatmap_final",width=10,height=6)

Heatmap(dat1, name='Global Fst', col = colorRamp2(c(0, 0.2), c("white",  "red")),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        cluster_rows = FALSE, cluster_columns = FALSE,
        rect_gp = gpar(type = "none"), 
        cell_fun = function(j, i, x, y, w, h, col) {
          if(i < j) {
            grid.rect(x, y, w, h, gp = gpar(fill = col))
          } else if(j == i) {
            #grid.text(labels[i], x, y)
          } else {
            grid.text(sprintf("%.2f", dat1[j, i]), x, y, gp = gpar(fontsize = 10))
          }
          grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "#E3E8F0"))
        })

graph2ppt(file="heatmap2",width=10,height=6)

############# make heatmap ############# 

set.seed(123)

mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

# permute the rows and columns 排列行和列
mat = mat[sample(nrow(mat), nrow(mat)), sample(ncol(mat), ncol(mat))]

rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)
Heatmap(mat, name = "foo", show_column_dend = FALSE)


p_noinbred <- c("#FF7F00","#FDBF6F","#fda56f","#871a1a","#E31A1c","#FB9A99","#8a584a" ,"#33A02C","#B2DF8A","#86a395","#1F78B4","#A6CEE3","#A6CEE3")
# Code below is way more complicated than it needs to be.  I'll share these files.
popmap$Type <- factor(popmap$Type, levels=c("Inbred","Domesticated","Wild"))
popmap$POP <- factor(popmap$POP, levels=c("HI","SM","UMFS","NEH","NG", "HG","HC","CS","DEBY" ,"CLP","HC-VA","LOLA","CL","SL","OBOYS","LM"))
