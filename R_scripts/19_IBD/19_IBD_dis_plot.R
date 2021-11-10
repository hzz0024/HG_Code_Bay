# https://boulter.com/gps/distance/?from=39.2935%09-75.3435&to=39.248%09-75.248&units=k
library(ggplot2)
library(cowplot)
IBD <- read.csv("IBD_dis_plot.csv", header = TRUE)
IBD$Fst

ggplot(data = IBD, aes(x = Distance, y = Fst)) +
  geom_point(shape=20, fill="red", color="black", size=5)+
  theme(axis.title.y=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  scale_y_continuous(name="Fst/(1-Fst)") +
  scale_x_continuous(name="Distance(km)", limits=c(0, 32))+
  geom_smooth(method = "lm", se = FALSE)
  
