setwd("/Volumes/cornell/CVreseq_plot/C_score/FigureYa83enrichment")
################ C_score test ################ 
d1<-read.csv("shared_outlier_noinvers_1K.csv", row.names = 1, header= TRUE)
#d1<-read.csv("C_score.csv", row.names = 1, header= TRUE)
d <- round((d1),2)
as.matrix(d)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(d){
  d[upper.tri(d)] <- NA
  return(d)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(d){
  d[lower.tri(d)]<- NA
  return(d)
}
upper_tri <- get_upper_tri(d)
upper_tri
lower_tri <- get_lower_tri(d)
lower_tri
rotate1 = t(apply(upper_tri, 2, rev))
rotate1
rotate2 = t(apply(rotate1, 2, rev))
rotate2
rotate = t(apply(rotate2, 2, rev))
rotate
# Melt the correlation matrix
library(reshape2)
melted_d <- melt(rotate, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_d, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "#FCC0C5", mid = "grey80", 
                       midpoint = 0, limit = c(0,60), space = "Lab", 
                       name= "Shared outliers") + #name=bquote("C"[hyper])) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

library(export)
graph2ppt(file="shared_outlier_noinvers_1K",width=10,height=6)

############# p-values ###############

d1<-read.csv("p-value.csv", row.names = 1, header= TRUE)
d <- round((d1),3)
as.matrix(d)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(d){
  d[upper.tri(d)] <- NA
  return(d)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(d){
  d[lower.tri(d)]<- NA
  return(d)
}
upper_tri <- get_upper_tri(d)
upper_tri
lower_tri <- get_lower_tri(d)
lower_tri
rotate1 = t(apply(upper_tri, 2, rev))
rotate1
rotate2 = t(apply(rotate1, 2, rev))
rotate2
rotate = t(apply(rotate2, 2, rev))
rotate
# Melt the correlation matrix
library(reshape2)
melted_d <- melt(rotate, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_d, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#E9EAE0", high = "#868B8E", mid = "#ADA7A7", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="p-value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

library(export)
graph2ppt(file="pvalue",width=10,height=6)

