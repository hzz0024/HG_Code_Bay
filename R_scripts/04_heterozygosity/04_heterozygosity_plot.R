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
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/04_heterozygosity")
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
library(wesanderson)
library(cowplot) 
 # OLD color values=c("#1BA3C6", "#33A65C", "#F8B620", "#E03426", "#EB73B3", "#AEC7E8", "#FF7F0E", "#9EDAE5", "#FFBB78")
cbPalette <- c("#A71B4B", "#E97302", "#EAC728", "#0BC0B3", "#4461A8", "#F5191C", "#7A7A7A")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/04_heterozygosity")
file = "DelBay_n_1478_het_summary.csv"
df <- read.delim(file, header = TRUE, sep=',')
df <- df[which(df$Heterozygosity > 1e-6),]
# reorder the populations # https://www.r-graph-gallery.com/22-order-boxplot-labels-by-names.html
#df$Pop <- factor(df$Pop , levels=c("HC", "ARN", "COH", "SR", "NB", "19 Surv.","19 Ref","20 Surv.","20 Ref"))
df$Pop <- factor(df$Pop , levels=c("18_HC", "18_ARN", "18_COH", "18_SR", "18_NB", 
                                   "19_HC", "19_ARN", "19_COH", "19_SR", "19_NB",
                                   "21_HC", "21_ARN", "21_COH", "21_SR", "21_BS", "21_BEN", "21_NAN", "21_NB", 
                                   "21_spat_HC", "21_spat_ARN", "21_spat_Nan", "19_Sur", "19_Ref", "20_Sur", "20_Ref", "A_20", "B_20", "C_20" ))

cbPalette <- wes_palette("Zissou1", 28, type = "continuous")
order <- c("18_HC", "18_ARN", "18_COH", "18_SR", "18_NB", 
           "19_HC", "19_ARN", "19_COH", "19_SR", "19_NB",
           "21_HC", "21_ARN", "21_COH", "21_SR", "21_BS", "21_BEN", "21_NAN", "21_NB", 
           "21_spat_HC", "21_spat_ARN", "21_spat_Nan", "19_Sur", "19_Ref", "20_Sur", "20_Ref", "A_20", "B_20", "C_20" )

p1 <- ggplot(df, aes(x=Pop, y=Heterozygosity, fill=as.factor(Pop))) +
  geom_boxplot(linetype = "dashed", outlier.shape = NA) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  #geom_violin(trim=FALSE, alpha = 0.8, width=1.1)+
  #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="black")+
  theme(legend.position="right")+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  #scale_colour_manual(values=cbPalette, breaks=order)+
  scale_fill_manual(values=cbPalette, breaks=order)+
  theme_classic()+ 
  guides(fill=guide_legend(title="Population"))+
  theme(text=element_text(family="Times New Roman", face="bold", size=14, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))

p1

graph2ppt(file="04_het_DelBay.pptx", width=14, height=10)

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
