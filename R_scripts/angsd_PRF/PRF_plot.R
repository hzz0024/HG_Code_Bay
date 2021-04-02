
####### plot the PRF accuracy for top SNPs #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd("~/Documents/Ryan_workplace/DelBay_adult/14_PRF")
file = "Del19_3006.txt"
df <- read.delim(file, header = TRUE, sep='\t')
as.factor(df$Top_SNPs)
ggplot(df, aes(x=Top_SNPs, y=Accuracy, shape=Top_SNPs)) +
  geom_point()

p1 <- ggplot(df, aes(x=Top_SNPs, y=Accuracy, fill=as.factor(Top_SNPs))) +
  ylim(0, 1) + # color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right") +
  geom_line(data = df, aes(x = Top_SNPs, y = Accuracy))
p1

graph2ppt(file="Del19_accuracy",width=18,height=9)


