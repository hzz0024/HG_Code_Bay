library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(export)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/26_salinity_index")
df<-read.delim("Correlation_salinityindex_mortality.txt", header = TRUE, sep='\t')

mytable <- df %>% 
  group_by(Year) %>%
  summarize(cor=cor(TotalMortality, Mean_salinity),
            p = cor.test(TotalMortality, Mean_salinity)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p1 <- df %>% 
  group_by(Year) %>% 
  ggplot(aes(x=TotalMortality, y=Mean_salinity, group=Year, color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)

p1

mytable <- df %>% 
  group_by(as.factor(Year)) %>%
  summarize(cor=cor(TotalMortality, CD5),
            p = cor.test(TotalMortality, CD5)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p2 <- df %>% 
  group_by(as.factor(Year)) %>% 
  ggplot(aes(x=TotalMortality, y=CD5, group=as.factor(Year), color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)

p2

mytable <- df %>% 
  group_by(as.factor(Year)) %>%
  summarize(cor=cor(TotalMortality, CD7),
            p = cor.test(TotalMortality, CD7)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p3 <- df %>% 
  group_by(as.factor(Year)) %>% 
  ggplot(aes(x=TotalMortality, y=CD7, group=as.factor(Year), color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)
p3

mytable <- df %>% 
  group_by(as.factor(Year)) %>%
  summarize(cor=cor(TotalMortality, CD9),
            p = cor.test(TotalMortality, CD9)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p4 <- df %>% 
  group_by(as.factor(Year)) %>% 
  ggplot(aes(x=TotalMortality, y=CD9, group=as.factor(Year), color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)

p4

mytable <- df %>% 
  group_by(as.factor(Year)) %>%
  summarize(cor=cor(TotalMortality, CD11),
            p = cor.test(TotalMortality, CD11)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p5 <- df %>% 
  group_by(as.factor(Year)) %>% 
  ggplot(aes(x=TotalMortality, y=CD11, group=as.factor(Year), color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)

p5

mytable <- df %>% 
  group_by(as.factor(Year)) %>%
  summarize(cor=cor(TotalMortality, MAX10),
            p = cor.test(TotalMortality, MAX10)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

p6 <- df %>% 
  group_by(as.factor(Year)) %>% 
  ggplot(aes(x=TotalMortality, y=MAX10, group=as.factor(Year), color=as.factor(Year))) +
  geom_point(size=6) +
  scale_color_manual(values=c("#ACD0C0", "#D09683", "#73605B"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  annotation_custom(tableGrob(mytable), xmin=45, xmax=60, ymin=0, ymax=-1)

p6

p1+p2+p3+p4+p5+p6+plot_layout(guides="collect")
graph2ppt(file="salinityindex_mortality_cor",width=18,height=12)

##################
# GLM regression #
##################

pop <- factor(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)), levels
              = c("Navarra", "Aragon", "Catalonia")) # Population
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9) # Wingspan
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2) # Body length
sex <- factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F"))
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6) # Number of ectoparasites
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57) # Color intensity
damage <- c(0,2,0,0,4,2,1,0,1) # Number of wings damaged
cbind(pop, sex, wing, body, mites, color, damage) # Print out data set
str(pop)
par(mfrow = c(1, 3), cex = 1.2)
colorM <- c("red", "red", "blue", "green", "green") # Pop color code males
colorF <- c("red", "blue", "blue", "green", "green") # Pop color code females
plot(body[sex == "M"], wing[sex == "M"], col =colorM, xlim = c(6.5, 9.5), ylim = c(10, 14),
     lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wingspan")
points(body[sex == "F"], wing[sex == "F"], col =colorF, pch = 16)
text(6.8, 13.8, "A", cex = 1.5)
plot(body[sex == "M"], mites[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 10),
     lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Parasite load")
points(body[sex == "F"], mites[sex == "F"], col = colorF, pch = 16)
text(6.8, 9.7, "B", cex = 1.5)
plot(body[sex == "M"], damage[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 4),
     lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Damaged wings")
points(body[sex == "F"], damage[sex == "F"], col = colorF, pch = 16)
text(6.8, 3.9, "C", cex = 1.5)



