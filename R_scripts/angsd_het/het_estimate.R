# the imputs are .ml files, which generated from the angsd
# #!/bin/sh
# source /local/workdir/hz269/DelBay_all_angsd/01_scripts/01_config.sh
# realSFS="/programs/angsd20191002/angsd/misc/realSFS"
# REGIONS="-rf chr_list.txt" #optional
# 
# for sample in `cat bam.list`; do
# SAMPLE_ID=`echo $sample | cut -d$'_' -f 1-2`
# echo "Sample is :$SAMPLE_ID "
# echo "#####################################"
# $angsd -i $sample -GL 1 -minQ 20 -minmapq 30 -anc $ANC -dosaf 1 -fold 1 $REGIONS -out $SAMPLE_ID
# $realSFS $SAMPLE_ID".saf.idx" > $SAMPLE_ID".ml"
# done



####### format data and produce the csv file #######
library(data.table)
setwd('~/Documents/Ryan_workplace/DelBay_adult/batch_effect/heterozygosity/het1')
setwd('~/Documents/Ryan_workplace/DelBay_adult/batch_effect/heterozygosity/het2')
files <- list.files(pattern="*.ml", full.names=TRUE, recursive=FALSE)
outs = sapply(files, function(x) {
  t <- read.delim(x, header = FALSE, sep=' ')
  out <- (t[2]/sum(t[1:2]))
  out = unlist(out)
  print(out)
})
write.table(outs, "./DelBay_minq33.csv", sep=",", quote=T, row.names=T, col.names=NA)

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd('~/Documents/Ryan_workplace/DelBay_adult/batch_effect/heterozygosity/plot/')
file = "DelBay_het_comp.csv"
file = "DelBay_het_chr.csv"
df <- read.delim(file, header = TRUE, sep=',')
ggplot(df, aes(x=Het_minq33, y=Het_minq20, shape=Pop)) +
  geom_point()

p1 <- ggplot(df, aes(x=Pop, y=Het_minq33, fill=as.factor(Batch))) +
  ylim(0.0056, 0.009) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="top")
p1
p2 <- ggplot(df, aes(x=Pop, y=Het_minq20, fill=as.factor(Batch))) + 
  ylim(0.0056, 0.009) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="top")
p2
grid.arrange(p1, p2, nrow = 1)
graph2ppt(file="comp",width=18,height=9)

############################### below are analysis across all sequenced bam files ##############################
####### format data and produce the csv file #######
library(data.table)
setwd("~/Documents/Ryan_workplace/DelBay_adult/04_heterozygosity/plot/")
files <- list.files(pattern="*.ml", full.names=TRUE, recursive=FALSE)
outs = sapply(files, function(x) {
  t <- read.delim(x, header = FALSE, sep=' ')
  out <- (t[2]/sum(t[1:2]))
  out = unlist(out)
  print(out)
})
write.table(outs, "./DelBay_minq20.csv", sep=",", quote=T, row.names=T, col.names=NA)

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd('~/Documents/Ryan_workplace/DelBay_adult/04_heterozygosity/plot/')
file = "DelBay_het_comp.csv"
file = "DelBay_het_chr.csv"
file = "DelBay_minq20_formal.csv"
df <- read.delim(file, header = TRUE, sep=',')
ggplot(df, aes(x=Het_minq20, y=Het_minq20, shape=Pop)) +
  geom_point()

p1 <- ggplot(df, aes(x=Pop, y=Het_minq20, fill=as.factor(Batch))) +
  #ylim(0.0056, 0.009) + color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p1
p2 <- ggplot(df, aes(x=Pop, y=Het_minq20, fill=as.factor(Batch))) + 
  ylim(0.0056, 0.009) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="top")
p2
grid.arrange(p1, p2, nrow = 1)
graph2ppt(file="all_441",width=18,height=9)
