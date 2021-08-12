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
setwd("~/Documents/HG/DelBay19_adult/04_heterozygosity/results")
files <- list.files(pattern="*.ml", full.names=TRUE, recursive=FALSE)
outs = sapply(files, function(x) {
  t <- read.delim(x, header = FALSE, sep=' ')
  out <- (t[2]/sum(t[1:2]))
  out = unlist(out)
  print(out)
})
write.table(outs, "./DelBay19_het_summary.csv", sep=",", quote=T, row.names=T, col.names=NA)

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)

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
write.table(outs, "./DelBay19_het_summary.csv", sep=",", quote=T, row.names=T, col.names=NA)

####### plot the het for different minq results #######
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
setwd("~/Documents/HG/DelBay19_adult/04_heterozygosity/plot")
file = "DelBay19_het_summary.csv"
df <- read.delim(file, header = TRUE, sep=',')

p1 <- ggplot(df, aes(x=Pop, y=Heterozygosity, fill=as.factor(Pop))) +
  #ylim(0.0056, 0.009) + color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="right")
p1

graph2ppt(file="04_het.pptx", width=10, height=6)

p2 <- ggplot(df, aes(x=Pop, y=Heterozygosity, fill=as.factor(Pop))) + 
  ylim(0.0, 0.009) + #color=as.factor(Batch), shape=as.factor(Batch)
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme(legend.position="top")
p2
grid.arrange(p1, p2, nrow = 1)
graph2ppt(file="all_441",width=18,height=9)
