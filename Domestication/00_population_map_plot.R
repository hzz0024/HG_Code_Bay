# ### Setup R Environment, Colors, and plotting set
library("dplyr")
library("ggplot2")
library(devtools)
library(R.utils)
library(gghighlight)
library(ggman)
library(ggtext)
library(patchwork)
library(plotrix)
library(plyr)
library(qqman)
library(qvalue)
library(reshape2)
library(tidyr)
library(zoo)
library(infer)
options(dplyr.summarise.inform = FALSE)
library(bigsnpr)
library("wesanderson")
library("directlabels")
library(OutFLANK)
library(pcadapt)
library(adegenet)
library(poppr)
library(vcfR)
library(maps)
library(mapdata)
library(fiftystater)
library("rgdal") 
library("plyr")
library(RColorBrewer)
library(matrixStats)
library(data.table)
library(ggdendro)
library(ggridges)
library(ggrepel)
library(legendMap)
library(export)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/00_pop_info/")

pops <- read.table("./Site_info_for_map.txt", header = TRUE)
# popind <- read.table("./00_pop_info/popmap.txt", header=TRUE)
# popmap <- pops %>% inner_join(popind,by="POP")
# popmap <- popmap[order(popmap$IND),]
# 
# # order population manually
# popmap$Type <- factor(popmap$Type, levels=c("Selected","Wild"))
# popmap$POP <- factor(popmap$POP, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW2", "NCW1", "MEH2", "NEH1", "NEH2", "UMFS", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
# 
# popmap.ni <- popmap # ni stands for non-inbreeding
#                 #   MEW1        MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
# p_noinbred <- c( "#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA")
# p1 <- p_noinbred
# 
# popmap.wild <- droplevels(subset(popmap.ni, POP != "MEH2" & POP !="NEH1" & POP != "NEH2" & POP != "UMFS" & POP != "DBX1" 
#                                  & POP != "DBX2" & POP != "DBX3" & POP != "UNC1" & POP != "UNC2"))
# p_wild <- c("#FF7F00","#FDBF6F","#E31A1c","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3")
# 
# gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
# gp.popmap$ind <- popmap.ni$IND
# gp.popmap$pop <- popmap.ni$POP
# 
# target <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW2", "NCW1", "MEH2", "NEH1", "NEH2", "UMFS", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2")
# gp.popmap <-gp.popmap %>% arrange(factor(pop, levels = target))

# Points to plot for each population and the population name
target <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2")

pop_coordinates$POP <- factor(pops$POP, levels = target)
pop_coordinates$Lat <- as.character(pop_coordinates$Lat)
pop_coordinates$Lat <- as.numeric(pop_coordinates$Lat)
pop_coordinates$Long <- as.character(pop_coordinates$Long)
pop_coordinates$Long <- as.numeric(pop_coordinates$Long)

        #   MEW1        MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
p1 <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd")

prefix.I <- "./Population_Genomics/Intermediate_Files/"
prefix.O <- "./Population_Genomics/Output/"
prefix.OE <- "./Population_Genomics/Output/Extra/"
prefix <- "./Population_Genomics/"

# pal1 <- wes_palette("Zissou1")
# pal2 <- wes_palette("GrandBudapest1")
# pal3 <- wes_palette("GrandBudapest2")
# pal4 <- wes_palette("Royal2")
# col_pal <- c(pal2[1],pal1[1], pal2[3],pal3[2],pal3[1],pal4[5],pal3[3],pal3[4],pal4[1],pal4[3])

# Map Creation

state <- map_data("state", boundary = FALSE)
w2hr <- map_data("world2Hires","USA")
states <- map_data("state")
# cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#999999","#F0E442" , "#CC79A7")
# pal <- wes_palette("Zissou1", 8, type = "continuous")

coastline = readOGR(dsn="./other_files/map_data", layer="us_medium_shoreline")
states= readOGR(dsn="./other_files/map_data", layer="state_bounds")
coast <- spTransform(coastline, CRS("+proj=longlat +ellps=clrk66"))
states <- spTransform(states, CRS("+proj=longlat +ellps=clrk66"))
coast@data$id = rownames(coastline@data)
states@data$id = rownames(states@data)
coastline.points = fortify(coastline, region="id")
states.points = fortify(states, region="id")

coast.df = join(coastline.points, coast@data, by="id")
states.df = join(states.points, coast@data, by="id")

filename <- paste(prefix.O,"PopulationMap.png", sep="")
png(filename=filename, type="cairo",units="px", width=2500, height=3500, res=200, bg="transparent")
p1 <- ggplot()+
  geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black" ,size=0.25, alpha=0.5) +
  geom_path(data = states.df, aes(x=long, y = lat, group = group), color="black",size=0.25, alpha=0.5)+
  geom_point(data=pop_coordinates, aes(x=Long,y=Lat, fill=POP, color=POP), size=8, shape=17,alpha=0.9)+
  scale_fill_manual(values=p1, name="Population")+
  scale_color_manual(values=p1, name="Population")+
  geom_label_repel(data=pop_coordinates, max.overlaps = Inf, aes(x=Long, y=Lat, label=POP, color=POP), fill="#f2f2f2", min.segment.length = 0,segment.colour = 'grey40',size=6, segment.size=0.5, box.padding=1,fontface="bold", alpha=1) + xlab("Longitude") + ylab("Latitude") + guides(fill="none", color="none")+
  xlab("Longitude") + 
  ylab("Latitude") +
  #coord_equal() +
  coord_map(xlim = c(min(pop_coordinates$Long), max(pop_coordinates$Long)),  ylim = c(min(pop_coordinates$Lat), max(pop_coordinates$Lat)), projection = "mercator")+
  theme_bw() +
  theme(text = element_text(size=24),plot.margin = margin(0, 0, 0, 0, "cm"))+
  scale_bar(lon = -72, lat = 34, distance_lon = 100, distance_lat = 10, distance_legend = 30, 
            dist_unit = "km", arrow_length = 50, arrow_distance = 60, arrow_north_size = 4)

graph2ppt(file="Map_wild",width=12,height=10)

dev.off()

# FST and Outlier Calculations


