# https://github.com/jorgeassis/marineDistances

## Source the main function 
source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")
## Read the landmass polygon
global.polygon <- "GSHHS_f_L1.shp"

## Run the function
setwd("~/Dropbox/CVreseq/ArcGIS")

## global.polygon: the path of the polygon
## file : the main file with the locations; should be text delimited
## file.strucutre: the main file structure: 1 to “Name Lon Lat” or 2 to “Name Lat Lon”
## file.header: define if the text file has a header with the column names (TRUE or FALSE)
## resolution: the resolution jjjjjjjjof the study area and the buffer to use around the sites. 
## buffer: the buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax). 
## export.file: file to export the results as a text delimited file (TRUE or FALSE)

contour(  global.polygon = "GSHHS_i_L1.shp" ,
          file = "Site_GPS.csv" , 
          file.sep = "," ,
          file.dec = "." ,
          file.strucutre = 2 , 
          file.header = TRUE ,
          resolution = 0.01,
          buffer = c(5,5,5,5) ,
          export.file = TRUE)

