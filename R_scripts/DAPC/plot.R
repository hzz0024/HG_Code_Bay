install.packages("export")
library(export)
library(ggplot2)
df2 <- read.csv("str.csv")
df2
# Line plot with multiple groups
p<-ggplot(df2, aes(x=K, y=N, group=Group)) +
  geom_line(aes(linetype=Group,color=Group), size = 1)+
  geom_point(aes(shape=Group, color=Group), size = 4)+
  scale_x_continuous(limits=c(1,8.1),breaks=c(1,2,4,6,8),labels = c("1","2","4","6","8"), expand = c(0, 0))+
  scale_y_continuous(limits=c(1,8), expand = c(0, 0))+
  scale_linetype_manual(values=c("twodash","longdash","solid","dotted"))+
  scale_shape_manual(values=c(15, 16, 17, 18))+
  scale_color_manual(values=c("#000000", "#474747", "#7f7f7f", "#c6c6c6")) + 
  theme_classic(base_size = 20)+
  xlab("Number of processor") + 
  ylab("Speed increase")+
  theme(legend.position = c(0.8, 0.8))+
  theme(legend.title = element_blank())
p

df1 <- read.csv("final_new.csv")
df1
# Line plot with multiple groups
p1<-ggplot(df1, aes(x=K, y=N, group=Group)) +
  geom_line(aes(linetype=Group,color=Group), size = 1)+
  geom_point(aes(shape=Group, color=Group), size = 4)+
  scale_x_continuous(limits=c(1,8.1),breaks=c(1,2,4,6,8),labels = c("1","2","4","6","8"), expand = c(0, 0))+
  scale_y_continuous(limits=c(1,8), expand = c(0, 0))+
  scale_linetype_manual(values=c("twodash","longdash","solid","dotted"))+
  scale_shape_manual(values=c(15, 16, 17, 18))+
  scale_color_manual(values=c("#000000", "#474747", "#7f7f7f", "#c6c6c6")) + 
  theme_classic(base_size = 20)+
  xlab("Number of processor") + 
  ylab("Speed increase")+
  theme(legend.position = c(0.8, 0.8))+
  theme(legend.title = element_blank())
p1

library(export)
graph2ppt(file="new",width=8,height=6)
