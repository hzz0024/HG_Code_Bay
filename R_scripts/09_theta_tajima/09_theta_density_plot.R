library(sm)
attach(mtcars)
cyl.f <- factor(cyl, levels= c(4,6,8),labels = c("4 cylinder", "6 cylinder", "8 cylinder"))
# plot densities
sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main="MPG Distribution by Car Cylinders")                

library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)

setwd("~/Documents/HG/DelBay19_adult/00_coverage/plot")

file = "coverage_het_plot.csv"
df <- read.delim(file, header = TRUE, sep=',')

p1 <- ggplot(df, aes(x=Pop1, y=Heterozygosity, fill=as.factor(Pop))) +
  #ylim(0.0056, 0.009) + color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p1
library(plyr)
sm.density.compare(df$Pop1, df$Mean_depth, xlab="Heterozygosity")
p <- ggplot(df, aes(x=Heterozygosity, color=Pop)) + 
  geom_density()
p
