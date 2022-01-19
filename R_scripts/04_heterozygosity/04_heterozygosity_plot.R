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
####################################################
####### format data and produce the csv file #######
####################################################
library(data.table)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/04_heterozygosity/Angsd_output")
files <- list.files(pattern="*.ml", full.names=TRUE, recursive=FALSE)
outs = sapply(files, function(x) {
  t <- read.delim(x, header = FALSE, sep=' ')
  out <- (t[2]/sum(t[1:2]))
  out = unlist(out)
  print(out)
})
write.table(outs, "./DelBay_het_summary.csv", sep=",", quote=T, row.names=T, col.names=NA)

####################################################
####### plot the het for different populations #####
####################################################
library("ggplot2")
library("Hmisc")
library("gtable")
library("gridExtra")
library(export)
library(cowplot) 
 # OLD color values=c("#1BA3C6", "#33A65C", "#F8B620", "#E03426", "#EB73B3", "#AEC7E8", "#FF7F0E", "#9EDAE5", "#FFBB78")
cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/04_heterozygosity/plot")
file = "DelBay19_het_summary.csv"
df <- read.delim(file, header = TRUE, sep=',')
# reorder the populations # https://www.r-graph-gallery.com/22-order-boxplot-labels-by-names.html
#df$Pop <- factor(df$Pop , levels=c("HC", "ARN", "COH", "SR", "NB", "19 Surv.","19 Ref","20 Surv.","20 Ref"))
df$Pop <- factor(df$Pop , levels=c("HC", "ARN", "COH", "SR", "NB", "Surv","Ref"))

p1 <- ggplot(df, aes(x=Pop, y=Heterozygosity, fill=as.factor(Pop))) +
  #geom_boxplot(fatten = NULL, alpha = 0.8)+
  geom_violin(trim=FALSE, alpha = 0.8, width=1)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="black")+
  theme(legend.position="right")+
  #scale_colour_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB", "Surv", "Ref"))+
  scale_fill_manual(values=cbPalette, breaks=c("HC", "ARN", "COH", "SR", "NB", "Surv", "Ref"))+
  cowplot::theme_cowplot()+
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        size = 0.5,
                                        linetype = "dashed"))

p1

graph2ppt(file="04_het_Del19.pptx", width=6, height=4)

####################################################
####### plot the het for different populations #####
####### after coverage equalization            #####
####################################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/04_heterozygosity/plot")

file = "DelBay_het_summary_downsampled.csv"
df <- read.delim(file, header = TRUE, sep='\t')
# reorder the populations # https://www.r-graph-gallery.com/22-order-boxplot-labels-by-names.html
df$Pop <- factor(df$Pop , levels=c("HC", "ARN", "COH", "SR", "NB", "19 Surv.","19 Ref","20 Surv.","20 Ref"))

p1 <- ggplot(df, aes(x=Pop, y=Heterozygosity, fill=as.factor(Pop))) +
  geom_boxplot(fatten = NULL, alpha = 0.8)+
  # plot the mean values
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, position=position_dodge(0.8))+
  theme_cowplot() + 
  theme(legend.position="right")+
  scale_fill_manual(values=c("#1BA3C6", "#33A65C", "#F8B620", "#E03426", "#EB73B3", "#AEC7E8", "#FF7F0E", "#9EDAE5", "#FFBB78")) 
   
p1
graph2ppt(file="04_het.pptx", width=10, height=6)
