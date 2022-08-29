############################
#  plot Bomaf distribution #
############################
library(ggpubr)
library(cowplot)
library(Rmisc)

#####################################
# Step 1: process the output files
#####################################
#setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/GEA_results")
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/results/")

format_outlier_output <- function(pop, j){
  #pop="HC_18"
  #j=1
  dat <- read.table(paste0("./outlier/", pop, "/", pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')

  output = data.frame(dat$chromo, dat$position, dat$CLR, dat$nSites, dat$SNP)
  colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  output = output[with(output, order(chromo, position)),]
  message(paste0("total number of snp is ", dim(output)[1]))
  #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
  return(output)
}

format_random_output <- function(pop, k, j){
  #pop="HC_18"
  #k=1
  #j=1
  if (file.exists(paste0("./random/", pop, "/", pop, "_", k, "_NC_03578", j ,".1.mafs.output.500bp.s1.txt"))) {
    dat <- read.table(paste0("./random/", pop, "/", pop, "_", k, "_NC_03578", j ,".1.mafs.output.500bp.s1.txt"), header = T)
    dat$chromo = paste0("NC_03578", j, '.1')
    dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
    dat = dat[with(dat, order(chromo, physPos)),]
    colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
    #random_list_name = "/Users/HG/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/18_HC_NB/18_HC_NB_random_2295"
    # random_list = paste0(random_list_name, "_", k, ".bed")
    # dat_ = read.delim(random_list, header = FALSE, sep='\t')
    # dat_$SNP = paste0(dat_$V1,'_',(dat_$V2+dat_$V3)/2)
  
    output = data.frame(dat$chromo, dat$position, dat$CLR, dat$nSites, dat$SNP)
    colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
    output = output[with(output, order(chromo, position)),]
    #message(paste0("total number of snp in ", k, "_", j, " is ", dim(output)[1]))
    #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
    return(output)
  } else {
    print("The file does not exist")
    return(data.frame())
  }
}

contrast <- function(pop, random_name){
  df_1=data.frame()
  for (j in seq(0,9)) {
    #pop="Sur_19"
    #outlier="18_HC_NB_shared_outlier.txt"
    df_1 = rbind(df_1,format_outlier_output(pop, j))
    colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  }
  message(paste0("total number of SNP in observed is ", dim(df_1)[1]))
  
  #df_1 <- df_1[which(df_1$CLR != 0), ]
  #df_1 <- df_1[which(df_1$nSites >= 10), ]
  
  df_random_mean = data.frame()
  for (i in seq(1,100)) {
    df_2 = data.frame()
    for (j in seq(0,9)){
      #random="/Users/HG/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/18_HC_NB/18_HC_NB_random_2295"
      df_2_ = format_random_output(pop, i, j)
      if(nrow(df_2_) != 0)
        df_2 = rbind(df_2, df_2_)
    }
    message(paste0("SNP counts in ", i, " is ", dim(df_2)[1]))
    #df_2 <- df_2[which(df_2$CLR != 0), ]
    #df_2 <- df_2[which(df_2$nSites >= 10), ]
    
    if(dim(df_1)[1] < dim(df_2)[1]){
      idx = seq(1:dim(df_2)[1])
      random = sort(sample(idx, dim(df_1)[1]))
      df_obs = df_1
      df_random = df_2[random,]
    } else {
      idx = seq(1:dim(df_1)[1])
      random = sort(sample(idx, dim(df_2)[1]))
      df_obs = df_1[random,]
      df_random = df_2
    }
    
    df_random_mean_ <- cbind(mean(df_obs$CLR, na.rm=TRUE), mean(df_random$CLR, na.rm=TRUE))
    df_random_mean <- rbind(df_random_mean, df_random_mean_)
  }
  
  df_random_mean_plot = as.data.frame(c(df_random_mean$V2))
  colnames(df_random_mean_plot) = c("random_mean")
  
  message("observed_mean: ", mean(df_random_mean$V1))
  message("random_mean_point: ", mean(df_random_mean_plot$random_mean))
  p <- sum(df_random_mean_plot$random_mean > mean(df_random_mean$V1))/dim(df_random_mean_plot)[1]
  message("p-value: ", p)
  return(df_random_mean_plot)
  #df_total = as.data.frame(c(df_random_mean$V1, df_random_mean$V2))
  #df_total$group = rep(c("outlier", "random"), c(dim(df_random_mean)[1], dim(df_random_mean)[1]))
  #colnames(df_total) = c("B0maf_score", "group")
}

HC_18 <- contrast("HC_18")
HC_18_observed_mean <- 29.5383515543699
HC_18_random_mean_point <- 27.9589291232822
HC_18_p <- 0
HC_19 <- contrast("HC_19")
HC_19_observed_mean <- 29.0325913203901
HC_19_random_mean_point <- 26.9239602803812
HC_19_p <- 0
HC_21 <- contrast("HC_21")
HC_21_observed_mean <- 27.1674500378662
HC_21_random_mean_point <- 25.0876191338405
HC_21_p <- 0
Sur_19 <- contrast("Sur_19")
Sur_19_observed_mean <- 24.0714660651283
Sur_19_random_mean_point <- 22.2598480243217
Sur_19_p <- 0
Sur_20 <- contrast("Sur_20")
Sur_20_observed_mean <- 27.482804884568
Sur_20_random_mean_point <- 25.5194265755036
Sur_20_p <- 0

NB_18 <- contrast("NB_18")
NB_18_observed_mean <- 21.6660260453865
NB_18_random_mean_point <- 23.3325268987198
NB_18_p <- 1
NB_19 <- contrast("NB_19")
NB_19_observed_mean <- 25.1534453003861
NB_19_random_mean_point <- 27.0982391762529
NB_19_p <- 1
NB_21 <- contrast("NB_21")
NB_21_observed_mean <- 23.7497553079883
NB_21_random_mean_point <- 25.4971145987965
NB_21_p <- 1
Ref_19 <- contrast("Ref_19")
Ref_19_observed_mean <- 21.3015528074266
Ref_19_random_mean_point <- 22.1307422886883
Ref_19_p <- 1
Ref_20 <- contrast("Ref_20")
Ref_20_observed_mean <- 44.6800864263947
Ref_20_random_mean_point <- 44.4468554077361
Ref_20_p <- 0.3
  
  
df_total = rbind(HC_18, HC_19, HC_21, Sur_19, Sur_20)
df_total$group = rep(c("HC_18", "HC_19", "HC_21", "Sur_19", "Sur_20"), c(dim(HC_18)[1], dim(HC_19)[1], dim(HC_21)[1], dim(Sur_19)[1], dim(Sur_20)[1]))
colnames(df_total) = c("random_mean", "group")
df_total$group = factor(df_total$group, levels = c("Sur_20" , "Sur_19", "HC_21", "HC_19", "HC_18"))

df_random = rbind(NB_18, NB_19, NB_21, Ref_19, Ref_20)
df_random$group = rep(c("NB_18", "NB_19", "NB_21", "Ref_19", "Ref_20"), c(dim(NB_18)[1], dim(NB_19)[1], dim(NB_21)[1], dim(Ref_19)[1], dim(Ref_20)[1]))
colnames(df_random) = c("random_mean", "group")
df_random$group = factor(df_random$group, levels = c("Ref_20" , "Ref_19", "NB_21", "NB_19", "NB_18"))

########
# plot #
########
mycolor<-c("#A71B4B", "#E97302", "#4461A8", "#0BC0B3","#EAC728")

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

mytheme<-theme_bw()+
  theme(
    axis.text.x = element_text(colour = "black", size=14),
    axis.text.y = element_text(colour = "black", size=14),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))

#一半小提琴图的参数调整：
#position：位置调整，这里将其向右水平移动0.1；
#side：显示哪一侧， "I"代表左侧，"R"代表右侧，默认"I"；
#adjust：调整带宽，这里设为1.2使宽带略变平滑；
#trim：小提琴图尾部的数据修整，默认为"T",表示将尾部修整到数据范围；"F"表示不修剪尾部；
p3_obs<-ggplot(df_total, aes(x=group, 
                             y=random_mean,
                             fill=group,
                             color=group))+
  scale_color_manual(values=(mycolor))+scale_fill_manual(values=(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  geom_point(aes(x = as.numeric(group)-0.1, 
                 y = random_mean,
                 color = group), 
             position = position_jitter(width = 0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  labs(y=expression(paste(Mean~B[0*', MAF'])), x = NULL)+
  geom_segment(aes(x = 5.4, y = HC_18_observed_mean, xend = 5.1, yend = HC_18_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 4.4, y = HC_19_observed_mean, xend = 4.1, yend = HC_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 3.4, y = HC_21_observed_mean, xend = 3.1, yend = HC_21_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 2.4, y = Sur_19_observed_mean, xend = 2.1, yend = Sur_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 1.4, y = Sur_20_observed_mean, xend = 1.1, yend = Sur_20_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  scale_linetype_manual(
    expression(paste(Observed~mean~B[0*', MAF'])), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p3_obs

p3_random<-ggplot(df_random, aes(x=group, 
                             y=random_mean,
                             fill=group,
                             color=group))+
  scale_color_manual(values=(mycolor))+scale_fill_manual(values=(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  geom_point(aes(x = as.numeric(group)-0.1, 
                 y = random_mean,
                 color = group), 
             position = position_jitter(width = 0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  labs(y=expression(paste(Mean~B[0*', MAF'])), x = NULL)+
  geom_segment(aes(x = 5.4, y = NB_18_observed_mean, xend = 5.1, yend = NB_18_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 4.4, y = NB_19_observed_mean, xend = 4.1, yend = NB_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 3.4, y = NB_21_observed_mean, xend = 3.1, yend = NB_21_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 2.4, y = Ref_19_observed_mean, xend = 2.1, yend = Ref_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 1.4, y = Ref_20_observed_mean, xend = 1.1, yend = Ref_20_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  scale_linetype_manual(
    expression(paste(Observed~mean~B[0*', MAF'])), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p3_random













# no longer useful
#########################
# Step 2: start to plot #
#########################
remotes::install_github("pstraforelli/tidytests")
library(gtools)
library(reshape2)
library(LDheatmap)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(Rmisc)
library(export)

df_plot <- df_total[df_total$group == "random",]
observed_mean <- mean(df_total[df_total$group == "outlier",]$B0maf_score)
random_mean_point <- mean(df_plot$B0maf_score)

ggplot(df_plot, aes(x=B0maf_score, fill=NULL)) +
  geom_histogram(fill="white", color="grey", bins = 15)+
  geom_vline(aes(xintercept=random_mean_point), color="grey",
             linetype="dashed")+
  labs(x=expression(Balancing~selection~score~"("~B0MAF~")"), y = "Count")+
  theme_classic() +
  ylim(0, 25)+
  theme(text = element_text(size = 15)) +
  geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
  geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
  geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 

