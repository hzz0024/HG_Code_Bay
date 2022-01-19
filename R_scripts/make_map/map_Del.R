install.packages("ggplot2")
install.packages("Rcpp")
install.packages("plyr")
install.packages("stringi")
install.packages("maptools")
install.packages("sp")
install.packages("Cairo")
install.packages("RColorBrewer")
install.packages("ggmap")
install.packages("sf")
install.packages("mapview")
library(plyr)
library(ggplot2)
library(sp)
library(maptools)
library(Cairo)
library(RColorBrewer)
library(ggmap)
library(sf)
library(mapview)

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/Map")
lat <- c(39.2,39.5)
long <- c(-75.6, -75.2)
bbox <- make_bbox(long,lat,f=0.05)
b <- get_map(bbox,maptype="toner", source = 'stamen')
ggmap(b)
#"#0000FF", "#ACEEF3", "#FFE9E4", "#FFB067", "#FF7077"
df <- read.csv("GPS_sites.csv")
head(df)
ggmap(b) + theme(axis.text=element_text(size = 20), axis.title = element_text(size = 15,)) +
  geom_point(data = df, aes(lon,lat,color=factor(Operator)), size = 6) + 
  scale_color_manual(values = c("#E97302", "#EAC728", "#A71B4B", "#4461A8", "#0BC0B3")) +
  labs(x = "Longitude", y = "Latitude", color = "") + 
  theme(legend.position="Null") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  cowplot::theme_cowplot()
graph2ppt(file="map", width=6, height=4)
