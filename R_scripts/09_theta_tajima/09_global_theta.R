setwd("~/Documents/HG/DelBay19_adult/09_theta/global_pi")
####### format data and produce the csv file #######
library(data.table)

files <- list.files(pattern="*.pestPG", full.names=TRUE, recursive=FALSE)
pi_cal = sapply(files, function(x) {
  t <- read.delim(x, header = TRUE, sep='\t')
  t$pi <- t$tP/t$nSites
  t = t[complete.cases(t), ]
  a = strsplit(x, split = "/")[[1]][2]
  pop = strsplit(a,split = ".th")[[1]][1]
  dat = cbind(t$Chr,t$pi, pop)
  print(dat)
})

write.table(pi_cal$`./ARN.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_ARN.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./COH.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_COH.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./HC.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_HC.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./NB.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_NB.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./SR.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_SR.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./REF19.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_REF19.csv", sep=",", quote=F, row.names=T, col.names=NA)
write.table(pi_cal$`./CHR19.thetas.1000.window.idx.pestPG`, "./DelBay19_pi_CHR19.csv", sep=",", quote=F, row.names=T, col.names=NA)

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd("~/Documents/HG/DelBay19_adult/09_theta/global_pi")
file = "DelBay19_pi_summary.csv"
df <- read.delim(file, header = TRUE, sep=',')

p1 <- ggplot(df, aes(x=Pop, y=Pi, fill=as.factor(Pop))) +
  ylim(0.0, 0.1) +
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p1

graph2ppt(file="04_het.pptx", width=10, height=6)

p2 <- ggplot(df, aes(x=Pop, y=Pi, fill=as.factor(Pop))) + 
  ylim(0, 0.02) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p2
grid.arrange(p1, p2, nrow = 1)
graph2ppt(file="DelBay19_global_pi",width=18,height=9)


##############################################


setwd("~/Documents/HG/DelBay19_adult/09_theta/global_pi_by_chr/")
####### format data and produce the csv file #######
library(data.table)

files <- list.files(pattern="*.pestPG", full.names=TRUE, recursive=FALSE)
pi_cal = sapply(files, function(x) {
  t <- read.delim(x, header = TRUE, sep='\t')
  t$pi <- t$tP/t$nSites
  t = t[complete.cases(t), ]
  a = strsplit(x, split = "/")[[1]][2]
  pop = strsplit(a,split = ".th")[[1]][1]
  dat = cbind(t$Chr,t$pi, pop)
  print(dat)
})

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)

file = "global_pi_by_chr.csv"
df <- read.delim(file, header = TRUE, sep=',')

p1 <- ggplot(df, aes(x=Pop, y=Pi, fill=as.factor(Pop))) +
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p1

graph2ppt(file="04_het.pptx", width=10, height=6)

p2 <- ggplot(df, aes(x=Pop, y=Pi, fill=as.factor(Pop))) + 
  ylim(0, 0.02) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p2
grid.arrange(p1, p2, nrow = 1)
graph2ppt(file="DelBay19_global_pi",width=18,height=9)


