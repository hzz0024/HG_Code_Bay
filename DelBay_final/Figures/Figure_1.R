#############################
# part 1 the geographic map #
#############################
library(raster)
library(maps) 
library(mapdata)
library(maptools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(ggsn)
library(tidyverse)
library(here)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_1")
BC.shp <- readOGR(here(".", "s_22mr22", "s_22mr22.shp"))
### chose the lat/long extent 
Ncalvert <- extent(-75.65, -75.15, 39.15, 39.50)
### crop your shapefile polygons to the extent defined
# takes a moment to run (patience grasshopper)
BC.shp2 <- crop(BC.shp,Ncalvert)
### project and fortify (i.e. turn into a dataframe)
BC.df <- fortify(BC.shp2)

DB_site <- read_csv("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_1/Site_GPS.csv")
DB_site$Location <- factor(DB_site$Location, levels = c("HC","ARN","COH","SR","BS","BEN","NAN","NB"))

p1<-ggplot()+ theme_bw()+
  geom_polygon(data = BC.df, aes(x=long,y=lat,group=group),
               colour= "black", size=0.2, fill="white")+
  coord_cartesian(xlim = c(-75.65, -75.15), ylim=c(39.15, 39.50)) +
  #geom_point(data=EXPTsites, aes(x=long, y=lat, shape=otter), size=4, colour="blue", stroke=1.5)+  #add this to plot site locations
  #scale_shape_manual(values=c(21,24))+         #this makes different shapes for otter "yes" and otter "no" sites
  scalebar(x.min = -75.65, x.max = -75.15, y.min = 39.15, y.max = 39.50,location = "bottomright",
           dist = 5,transform = TRUE, dist_unit = "km", model = 'WGS84', st.size = 4, st.dist = 0.02)+
  north(data = BC.df, scale = 0.1, symbol = 3, anchor= c(x = -75.15, y = 39.49)) +
  labs(shape="Year", colour="Population")+
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=12), axis.text.x = element_text(size=12))+
  geom_point(data=DB_site, size=3, mapping = aes(x = Lon, y = Lat, col= Location, shape=as.factor(Year))) +
  scale_color_brewer(palette = "Set2")

#############################
# part 2 the motarlity data #
#############################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/00_info/Mortality")
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(ggspectra)
library(magrittr)
library(tidyverse)
library(reshape2)
library(lubridate)
library(gridExtra)
library(patchwork)

input = "Figure1_mortality.txt"
#input = "short_term_monitoring.txt"
dat = read.delim(input, header = TRUE, sep='\t')
dat$Month[dat$Month == 1] = '01'
dat$Month[dat$Month == 2] = '02'
dat$Month[dat$Month == 3] = '03'
dat$Month[dat$Month == 4] = '04'
dat$Month[dat$Month == 5] = '05'
dat$Month[dat$Month == 6] = '06'
dat$Month[dat$Month == 7] = '07'
dat$Month[dat$Month == 8] = '08'
dat$Month[dat$Month == 9] = '09'
dat$yearmon <- paste0(dat$Year, "_", dat$Month)
head(dat)
dim(dat)
dat$Date = mdy(dat$Date)
new = levels(factor(dat$yearmon))
new_order = new[order(as.numeric(gsub("(\\d+)_(\\d+)", "\\1\\2", new)))]
dat$yearmon <- factor(dat$yearmon, levels = new_order)
cbPalette <- RColorBrewer::brewer.pal(8, "Set2")
#c("#4461A8", "#0BC0B3", "#EAC728", "#E97302", "#A71B4B")
dat$Population <- factor(dat$Population, levels = c("HC","ARN","COH","SR","BS","BEN","NAN_","NB"))

dat$Date <- as.POSIXct(dat$Date, format = "%Y-%m-%d")

df_plot <- data.frame(dat$Year, dat$Date, dat$Population, dat$Salinity, dat$TotalMortality)
colnames(df_plot) <- c("Year", "Date", "Population", "Salinity", "TotalMortality")
df_melt <- melt(df_plot, id = c("Year","Date", "Population"))
# you first create a named character vector for facet_grid https://stackoverflow.com/questions/3472980/how-to-change-facet-labels
variable_names <- c(
  `Salinity` = "Salinity",
  `TotalMortality` = "Mortality (%)"
)

p2<-ggplot(df_melt, aes(x=Date, y=value, group=Population, color = Population)) +
  geom_point( alpha=0.6, size=1) +
  geom_line(alpha=0.5, size=0.8) +
  #geom_smooth(method = loess, aes(fill = Population), alpha=0.4) +
  scale_shape_manual(rep(c(15),8), breaks=c("HC","ARN","COH","SR","BS","BEN","NAN_","NB")) + 
  scale_colour_manual(values=cbPalette, breaks=c("HC","ARN","COH","SR","BS","BEN","NAN_","NB"))+
  scale_fill_manual(values=cbPalette, breaks=c("HC","ARN","COH","SR","BS","BEN","NAN_","NB"))+
  theme_classic()+
  scale_x_datetime(date_labels="%Y", date_breaks = "1 year", expand=c(0,0)) +
  #stat_valleys(span=24, ignore_threshold = 10, color="blue") +
  #stat_valleys(geom="text", span=24, ignore_threshold = 10, x.label.fmt = "%Y", color="blue", angle=90, hjust=1.1) +
  theme(axis.text.y= element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_color_brewer(palette = "Set2")+
  facet_grid(variable ~ ., scales="free", labeller = as_labeller(variable_names))

p2

##############
# part 3 PCA #
##############

library("ggplot2")
library("cowplot")
library("ggrepel")
library("scales")
library("MASS")
library(export)
library("viridis")
library(wesanderson)

##########################################
# PCA for LD-pruned SNPs -- 1342 samples #
##########################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/10_prune_pca_mds")
source("individual_pca_functions.R")

PCA <- function(cov_matrix, ind_label, pop_label, x_axis, y_axis, show.point=T, show.label=T, show.ellipse=T, show.line=T, alpha=0, index_exclude=vector())
{
  ## This function takes a covariance matrix and performs PCA. 
  # cov_matrix: a square covariance matrix generated by most pca softwares
  # ind_label: a vector in the same order and length as cov_matrix; it contains the individual labels of the individuals represented in the covariance matrix
  # pop_label: a vector in the same order and length as cov_matrix; it contains the population labels of the individuals represented in the covariance matrix
  # x_axis: an integer that determines which principal component to plot on the x axis
  # y_axis: an integer that determines which principal component to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  # index_exclude: the indices of individuals to exclude from the analysis
  index_include <- setdiff(seq_along(ind_label), index_exclude)
  #m <- as.matrix(cov_matrix)
  m <- as.matrix(read.table(cov_matrix))
  m[is.na(m)]<- median(m, na.rm = T)
  m<-m[index_include, index_include] ## Remove 4SJH, four 3Ps individuals, and contaminated ones
  e <- eigen(m)
  e_value<-e$values
  x_variance<-e_value[x_axis]/sum(e_value)*100
  y_variance<-e_value[y_axis]/sum(e_value)*100
  e <- as.data.frame(e$vectors)
  e <- cbind(ind_label[index_include], pop_label[index_include], e) ## with the above individuals removed
  #colnames(e)[3:331]<-paste0("PC",1:329)
  colnames(e)[3:(dim(e)[1])]<-paste0("PC",1:(dim(e)[1]-2)) ## with the above individuals removed
  colnames(e)[1:2]<-c("individual", "population")
  assign("pca_table", e, .GlobalEnv)
  cbPalette <- wes_palette("Zissou1", 25, type = "continuous")
  pop_list <- c("18_HC","18_ARN","18_COH","18_SR","18_NB",
                "19_HC","19_ARN","19_COH","19_SR","19_NB",
                "21_HC","21_ARN","21_COH","21_SR","21_BS","21_BEN","21_NAN","21_NB",
                "19_Sur","19_Ref","20_Sur","20_Ref","20A","20B","20C")
  PCA_plot<-ggplot(data=e[,],aes(x=e[,x_axis+2], y=e[,y_axis+2], color=population,label=population, shape=population)) + 
    scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15), breaks=pop_list) +
    scale_colour_manual(values=cbPalette, breaks=pop_list)+
    geom_enterotype(alpha=alpha, show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line) +
    #geom_text_repel()+
    theme_bw()+
    #theme_cowplot() +
    #geom_point(size = 6, alpha = 0.5)+ # @HG change the point size
    theme(
      axis.text.y= element_text(size=12), 
      axis.text.x = element_text(size=12)
    ) +
    xlab(paste0("PC", x_axis, "(",round(x_variance,2),"%)")) +
    ylab(paste0("PC", y_axis ,"(",round(y_variance,2),"%)")) 
  # function to plot the labels note for population with one individual it will report warning message
  # PCA_plot<-ggplot(data=e[,],aes(x=e[,x_axis+2], y=e[,y_axis+2], color=population,label=individual, shape=population)) + 
  #   geom_enterotype(alpha=alpha, show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line) +
  #   scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  #   geom_text_repel()+
  #   theme_cowplot() +
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.text.y = element_blank()
  #   ) +
  #   xlab(paste0("PC", x_axis, "(",round(x_variance,2),"%)")) +
  #   ylab(paste0("PC", y_axis ,"(",round(y_variance,2),"%)")) 
  print(PCA_plot)
}

# control groups
file = 'All_1342_info.txt' 
all_label = read.delim(file, header = TRUE, sep='\t')
all_label$pop <- factor(all_label$pop, levels = c("18_HC","18_ARN","18_COH","18_SR","18_NB",
                                                  "19_HC","19_ARN","19_COH","19_SR","19_NB",
                                                  "21_HC","21_ARN","21_COH","21_SR","21_BS","21_BEN","21_NAN","21_NB",
                                                  "19_Sur","19_Ref","20_Sur","20_Ref","20A","20B","20C"))
# cbPalette <- RColorBrewer::brewer.pal(25, "Set2")
cbPalette <- wes_palette("Zissou1", 25, type = "continuous")

p3 <- PCA(cov_matrix = "All_pca_1342_all_minq20_minmq30_CV30_masked_chr9.covMat",
          ind_label = all_label$ind,
          pop_label = all_label$pop,
          x_axis = 1,
          y_axis = 2,
          show.point = T,
          show.label = F,
          show.ellipse = T,
          show.line = F,
          alpha = 0.01,
)

p3+ guides(shape = guide_legend(override.aes = list(size = 4)))+
  theme(legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7))

layout <-" 
AAAAAAAABBBBBBB
AAAAAAAABBBBBBB
AAAAAAAABBBBBBB
DDDDDDDDCCCCCCC
DDDDDDDDCCCCCCC
"
p1 + p2 + p3 + plot_layout(design=layout, guides='keep')

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_1")
graph2ppt(file="Figure_1.pptx", width=18, height=12)
