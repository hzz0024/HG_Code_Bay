library(gtools)
library(reshape2)
library(LDheatmap)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(Rmisc)
library(export)
####################################
##########  plot r2 values  ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/random_SNPs_nomatch")

format <- function (pname1, pname2, i){
  #pname1 = "./outlier_output/SR_HC/HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./random_output/HC.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < 500),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
  }
  df_total = rbind(dat1, dat2)
  df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  return(df_total)
}

chr1 <- format("./outlier_output/SR_HC/HC.chr1.ngsld.output", "./random_output/SR_HC/HC.chr1.ngsld.output", 1)
chr2 <- format("./outlier_output/SR_HC/HC.chr2.ngsld.output", "./random_output/SR_HC/HC.chr2.ngsld.output", 2)
chr3 <- format("./outlier_output/SR_HC/HC.chr3.ngsld.output", "./random_output/SR_HC/HC.chr3.ngsld.output", 3)
chr4 <- format("./outlier_output/SR_HC/HC.chr4.ngsld.output", "./random_output/SR_HC/HC.chr4.ngsld.output", 4)
chr5 <- format("./outlier_output/SR_HC/HC.chr5.ngsld.output", "./random_output/SR_HC/HC.chr5.ngsld.output", 5)
chr6 <- format("./outlier_output/SR_HC/HC.chr6.ngsld.output", "./random_output/SR_HC/HC.chr6.ngsld.output", 6)
chr7 <- format("./outlier_output/SR_HC/HC.chr7.ngsld.output", "./random_output/SR_HC/HC.chr7.ngsld.output", 7)
chr8 <- format("./outlier_output/SR_HC/HC.chr8.ngsld.output", "./random_output/SR_HC/HC.chr8.ngsld.output", 8)
chr9 <- format("./outlier_output/SR_HC/HC.chr9.ngsld.output", "./random_output/SR_HC/HC.chr9.ngsld.output", 9)
chr10 <- format("./outlier_output/SR_HC/HC.chr10.ngsld.output", "./random_output/SR_HC/HC.chr10.ngsld.output", 10)


df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

# Pairwise comparisons
pwc <- df %>% 
  pairwise_t_test(
    V7 ~ group, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

f1 <- function(x){
  log10(mean(10^x))
}
# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "group", step.increase = 1)
df1 = df[df$V7!=0,]
df2 <- summarySE(df1, measurevar="V7", groupvars=c("group"))
#df1$log = -log(df1$V7, 10)
ggplot(df1, 
       aes(x = group, y = V7)) + 
  stat_summary(fun = mean, geom = "pointrange") +
  stat_summary(fun = mean, geom = "line")


ggplot(df2, aes(x = group, y = V7)) + 
  geom_errorbar(width=.1, aes(ymin=V7-ci, ymax=V7+ci)) +
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1))+ 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6, vjust = -0.5)+
  guides(fill = "none") +
  #ylim(-10,10) +
  labs(x = NULL, y = expression(r^2)) +
  theme(text = element_text(size = 15)) +
  #coord_trans(y='log10') +
  coord_cartesian(ylim=c(0.15, 0.2))
  #coord_cartesian(ylim=c(0.15, 0.2))

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(df2, aes(x=group, y=V7, colour=group, group=group)) + 
  geom_errorbar(aes(ymin=V7-se, ymax=V7+se), colour="black", width=.1, position=pd) +
  coord_cartesian(ylim=c(0.17, 0.18))+
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

graph2ppt(file="HC-NB.pptx", width=4, height=3) 

############## 2019 challenge ###########

format <- function (pname1, pname2, i){
  #pname1 = "./outlier_output/SR_HC/HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./random_output/HC.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < 500),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
  }
  df_total = rbind(dat1, dat2)
  df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  return(df_total)
}

chr1 <- format("./outlier_output/REF19_CHR19/CHR19.chr1.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr1.ngsld.output", 1)
chr2 <- format("./outlier_output/REF19_CHR19/CHR19.chr2.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr2.ngsld.output", 2)
chr3 <- format("./outlier_output/REF19_CHR19/CHR19.chr3.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr3.ngsld.output", 3)
chr4 <- format("./outlier_output/REF19_CHR19/CHR19.chr4.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr4.ngsld.output", 4)
chr5 <- format("./outlier_output/REF19_CHR19/CHR19.chr5.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr5.ngsld.output", 5)
chr6 <- format("./outlier_output/REF19_CHR19/CHR19.chr6.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr6.ngsld.output", 6)
chr7 <- format("./outlier_output/REF19_CHR19/CHR19.chr7.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr7.ngsld.output", 7)
chr8 <- format("./outlier_output/REF19_CHR19/CHR19.chr8.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr8.ngsld.output", 8)
chr9 <- format("./outlier_output/REF19_CHR19/CHR19.chr9.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr9.ngsld.output", 9)
chr10 <- format("./outlier_output/REF19_CHR19/CHR19.chr10.ngsld.output", "./random_output/REF19_CHR19/CHR19.chr10.ngsld.output", 10)

df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

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
  coord_cartesian(ylim=c(0.18, 0.19))+
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

graph2ppt(file="2019_challenge.pptx", width=4, height=3) 

############## NB-HC ###########

format <- function (pname1, pname2){
  #pname1 = "./outlier_output/SR_HC/HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./random_output/HC.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < 500),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
  }
  df_total = rbind(dat1, dat2)
  df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  return(df_total)
}

df <- format("random_output/NB_HC/HC.2.ngsld.output", "./random_output/NB_HC/HC.1.ngsld.output")
df <- format("outlier_output/NB_HC/HC.outlier.ngsld.output", "./random_output/NB_HC/HC.1.ngsld.output")
chr1 <- format("./outlier_output/NB_HC/HC.chr1.ngsld.output", "./random_output/NB_HC/HC.chr1.ngsld.output", 1)
chr2 <- format("./outlier_output/NB_HC/HC.chr2.ngsld.output", "./random_output/NB_HC/HC.chr2.ngsld.output", 2)
chr3 <- format("./outlier_output/NB_HC/HC.chr3.ngsld.output", "./random_output/NB_HC/HC.chr3.ngsld.output", 3)
chr4 <- format("./outlier_output/NB_HC/HC.chr4.ngsld.output", "./random_output/NB_HC/HC.chr4.ngsld.output", 4)
chr5 <- format("./outlier_output/NB_HC/HC.chr5.ngsld.output", "./random_output/NB_HC/HC.chr5.ngsld.output", 5)
chr6 <- format("./outlier_output/NB_HC/HC.chr6.ngsld.output", "./random_output/NB_HC/HC.chr6.ngsld.output", 6)
chr7 <- format("./outlier_output/NB_HC/HC.chr7.ngsld.output", "./random_output/NB_HC/HC.chr7.ngsld.output", 7)
chr8 <- format("./outlier_output/NB_HC/HC.chr8.ngsld.output", "./random_output/NB_HC/HC.chr8.ngsld.output", 8)
chr9 <- format("./outlier_output/NB_HC/HC.chr9.ngsld.output", "./random_output/NB_HC/HC.chr9.ngsld.output", 9)
chr10 <- format("./outlier_output/NB_HC/HC.chr10.ngsld.output", "./random_output/NB_HC/HC.chr10.ngsld.output", 10)

df <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

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
  coord_cartesian(ylim=c(0.17, 0.19))+
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

graph2ppt(file="HC-NB.pptx", width=4, height=3) 


# ARN_COH
df <- format("outlier_output/ARN_COH/ARN.ngsld.output", "./random_output/ARN_COH/ARN.1.ngsld.output")
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
  coord_cartesian(ylim=c(0.165, 0.18))+
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

graph2ppt(file="ARN_COH.pptx", width=4, height=3) 

# REF19_SR
df <- format("./outlier_output/REF19_SR/REF19.ngsld.output", "./random_output/REF19_SR/REF19.1.ngsld.output")
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
  coord_cartesian(ylim=c(0.17, 0.19))+
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

graph2ppt(file="REF19_SR.pptx", width=4, height=3) 


####################################
##########  plot LD block   ########
####################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD")
TMP_FILE = "NB.ngsld.output"

df <- read.table(TMP_FILE, header=FALSE, stringsAsFactors=FALSE)
colnames(df)=c('snp1', 'snp2', 'dis', 'r1', 'D1', 'D2', 'r')
r <- df[complete.cases(df), ]
id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
posStart <- head(id,1)
posEnd <- tail(id,1)
r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))

SNPs = read.delim("chr5_ADF1.txt", header = FALSE, sep='\t')

for (ld in c("r")) {
  m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
  rownames(m) <- colnames(m)
  m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
  id <- rownames(m)
  dist <- as.numeric(sub(".*:","",id))
  
  # Save plot
  tiff(paste("LD_blocks", ld,"jpg", sep="."), units="in", width=6, height=6, res=300)
  #pdf(paste("LD_blocks", ld,"pdf", sep="."), width=10, height=10)
  LDheatmap(m, genetic.distances=dist, geneMapLabelX=0.75, geneMapLabelY=0.25, title = NULL, color="blueToRed", LDmeasure=ld, SNP.name = SNPs$V1)
  require(grid)
  grid.edit("symbols", pch = 20, gp = gpar(cex = 1, col = "red")) # change the symbol size
  grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=0.2, col = "red"))  # change the SNP label size
  dev.off()
}

#########################################
##########  plot LD distribution ########
#########################################
library(ggpubr)
library(cowplot)
format <- function (pname1, pname2){
  #pname1 = "./outlier_output/SR_HC/HC.chr1.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./random_output/HC.chr1.ngsld.output"
  dat2 = read.delim(pname2, header = FALSE, sep='\t')
  dat2 = dat2[,c(3,7)]
  dat2 = na.omit(dat2)
  dat2 = dat2[which(dat2$V3 < 500),]
  if(dim(dat1)[1] < dim(dat2)[1]){
    idx = seq(1:dim(dat2)[1])
    random = sort(sample(idx, dim(dat1)[1]))
    dat2 = dat2[random,]
  } else {
    idx = seq(1:dim(dat1)[1])
    random = sort(sample(idx, dim(dat2)[1]))
    dat1 = dat1[random,]
  }
  df_total = rbind(dat1, dat2)
  df_total$group = rep(c("outlier", "random"), c(dim(dat1)[1], dim(dat2)[1]))
  return(df_total)
}

random_mean <- c()
for(i in seq(1,100)){
  df <- format("outlier_output/NB_HC/HC.outlier.ngsld.output", paste0("./random_output/NB_HC/HC.", i, ".ngsld.output"))
  df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))
  message(i)
  print(df2$V7[2])
  random_mean = c(random_mean, df2$V7[2])
}

df<-format("./outlier_output/NB_HC/HC.outlier.ngsld.output", "./random_output/NB_HC/HC.1.ngsld.output")
observed_mean  <- summarySE(df, measurevar="V7", groupvars=c("group"))$V7[1]
random_mean <- as.data.frame(random_mean)
random_mean_point <- mean(random_mean$random_mean)

hist_p <- ggplot(random_mean, aes(x=random_mean, fill=NULL)) +
  geom_histogram(fill="white", color="grey", bins = 15)+
  geom_vline(aes(xintercept=mean(random_mean)), color="grey",
             linetype="dashed")+
  labs(x=expression(paste(Mean~Linkage~disequilibrium~(r^2))), y = "Count")+
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
  geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
  geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 
  
hist_p

graph2ppt(file="NB_HC_r2_comp", width=4, height=3) 

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in")
random_mean <- c()
for(i in seq(1,100)){
  df <- format("./outlier_output/CHR19.chr5.ngsld.output", paste0("./random_output/CHR19.", i, ".ngsld.output"))
  df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))
  message(i)
  print(df2$V7[2])
  random_mean = c(random_mean, df2$V7[2])
}

df<-format("./outlier_output/CHR19.chr5.ngsld.output", "./random_output/CHR19.1.ngsld.output")
observed_mean  <- summarySE(df, measurevar="V7", groupvars=c("group"))$V7[1]
random_mean <- as.data.frame(random_mean)
random_mean_point <- mean(random_mean$random_mean)

sum(random_mean$random_mean > observed_mean)

hist_p <- ggplot(random_mean, aes(x=random_mean, fill=NULL)) +
  geom_histogram(fill="white", color="grey", bins = 15)+
  geom_vline(aes(xintercept=mean(random_mean)), color="grey",
             linetype="dashed")+
  labs(x=expression(paste(Mean~Linkage~disequilibrium~(r^2))), y = "Count")+
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
  geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
  geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 

hist_p

graph2ppt(file="CHR19_REF19_r2_comp", width=4, height=3) 
