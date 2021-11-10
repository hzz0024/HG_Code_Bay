library(gtools)
library(reshape2)
library(LDheatmap)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

############################################
##########  format for random SNPs  ########
############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in/format")
format_bed <- function(pname, distance){
  #pname = "QTL_shared_loci"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  bed_list <- paste0(dat$V1, "\t" , dat$V2-distance, "\t", dat$V2+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("QTL_shared_loci", 250)
format_bed("QTL_random_loci", 250)

####################################
##########  format for bed  ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in/format")
format_bed <- function(pname, distance){
  #pname = "QTL_shared_loci"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  bed_list <- paste0(dat$V1, "\t" , dat$V2-distance, "\t", dat$V2+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("QTL_shared_loci", 250)
format_bed("QTL_random_loci", 250)
########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################

format_rf <- function(pname){
  dat = read.delim(pname, header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("QTL_shared_loci.bed.merged.txt")
format_rf("QTL_random_loci.bed.merged.txt")

####################################
##########  plot r2 values  ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in")
win = 500
format <- function (pname1, pname2, i){
  #pname1 = "./outlier_output/HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < win),]
  #pname2 = "./random_output/HC.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < win),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
    df_total = rbind(dat1, dat2)
    df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
    df_total = rbind(dat1, dat2)
    df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  }
  return(df_total)
}



#set.seed(2)
df <- format("./outlier_output/CHR19.chr5.ngsld.output", "./random_output/CHR19.chr5.ngsld.output", 1)

# Pairwise comparisons
pwc <- df %>% 
  pairwise_t_test(
    V7 ~ group, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "group", step.increase = 1)
ggboxplot(df, x = "group", y = "V7") +
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  ylim(0,1) +
  labs( x = NULL, y = expression(r^2))+
  theme(text = element_text(size = 15)) 

df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))
df2
ggplot(df2, aes(x=group, y=V7, colour=group, group=group)) + 
  geom_errorbar(aes(ymin=V7-se, ymax=V7+se), colour="black", width=.1, position=pd) +
  coord_cartesian(ylim=c(0.1, 0.45))+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1))+ 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6, hjust = 1.4, vjust = -0)+
  #geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  labs(x = NULL, y = expression(r^2)) +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("OJ", "VC"),
                   labels=c("Orange juice", "Ascorbic acid"),
                   l=40) +                    # Use darker colors, lightness=40
  theme_bw() +
  theme(text = element_text(size = 15)) 

graph2ppt(file="CHR19_REF19.pptx", width=4, height=3)


#########################
# B0maf plot for zoom in#
#########################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/zoom_in")
df_total = rbind(df_1, df_2)
df_total$group = rep(c("outlier", "random"), c(dim(df_1)[1], dim(df_2)[1]))
colnames(df_total)=c('chromo', 'position', 'B0maf_score', 'nSites', 'SNP', 'group')
df_total <- read.table("chr5_16552716.txt", header = T)
pwc <- df_total %>% 
  pairwise_t_test(
    B0maf_score ~ group, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

ggboxplot(df_total, x = "group", y = "B0maf_score") +
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  labs( x = NULL, y = expression(Balancing~selection~score~'('~B[0]~','~MAF~')'))+
  theme(text = element_text(size = 15)) 

df2 <- summarySE(df_total, measurevar="B0maf_score", groupvars=c("group"))
df2
ggplot(df2, aes(x=group, y=B0maf_score, colour=group, group=group)) + 
  geom_errorbar(aes(ymin=B0maf_score-se, ymax=B0maf_score+se), colour="black", width=.1, position=pd) +
  coord_cartesian(ylim=c(3, 12))+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1))+ 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6, hjust = 1.4, vjust = -0)+
  #geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  labs(x = NULL, y = expression(r^2)) +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("OJ", "VC"),
                   labels=c("Orange juice", "Ascorbic acid"),
                   l=40) +                    # Use darker colors, lightness=40
  theme_bw() +
  theme(text = element_text(size = 15)) 

graph2ppt(file="CHR19_REF19_zoomin_chr5.pptx", width=4, height=3)
