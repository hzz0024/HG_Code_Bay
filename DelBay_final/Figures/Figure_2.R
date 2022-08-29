
#####################
# plot for figure 2 #
#####################

# part1 sfs contrasts
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/20_balancing_selection/sfs/")

#function to normalize
norm <- function(x) x/sum(x)

spec_format <- function(sfs_name){
  sfs<-(scan(sfs_name))
  sfs_<-norm(sfs[-c(1,length(sfs))])
  barplot(sfs_) #plot variable sites 
  df <- data.frame(seq(1:length(sfs_)), length(sfs_)+1, sfs_)
  write.table(df, file = paste0("./spect_results/",sfs_name, ".nosub_spect.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

spec_format("NB_18.neutral.sfs")
spec_format("NB_19.neutral.sfs")
spec_format("NB_21.neutral.sfs")
spec_format("Ref_19.sfs")
spec_format("Ref_20.sfs")
spec_format("HC_18.neutral.sfs")
spec_format("HC_19.neutral.sfs")
spec_format("HC_21.neutral.sfs")
spec_format("Sur_19.sfs")
spec_format("Sur_20.sfs")

###################################
############ part 1 ###############
###################################
# https://r-graph-gallery.com/101_Manhattan_plot.html 
#BiocManager::install("karyoploteR")
###################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/12_shared_outlier_manhattan")
library(dplyr)
library(qqman)
library(karyoploteR)
library(patchwork)
library(ggplot2)

manhattan_plot <- function(pname, outlier_name, outlier_color, x_label){
  #pname="ps_18_SGS_HC_NB.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat$SNP = paste0(dat$V1 , "_", dat$V2)
  colnames(dat) = c("Chr", "Pos", "AF", "BF", "deltap", "ps", "signif", "SNP")
  # A list of SNP of interest is provided
  #outlier_name="18_HC_NB_shared_outlier.txt"
  outlier_file <- read.delim(outlier_name, header = T, sep='\t')
  snpsOfInterest <- outlier_file$id
  
  don <- dat %>% 
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len=max(Pos)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(dat, ., by=c("Chr"="Chr")) %>%
    # Add a cumulative position of each SNP
    arrange(Chr, Pos) %>%
    mutate(BPcum=Pos+tot) %>%
    mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) 
  
  axisdf = don %>% group_by(Chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  sig_data <- don %>% 
    subset(is_highlight == "yes")
  notsig_data <- don %>% 
    subset(is_highlight == "no") %>%
    group_by(Chr) %>% 
    sample_frac(0.01)
  
  plot_dat <- bind_rows(sig_data, notsig_data)
  
  ggplot(plot_dat, aes(x=BPcum, y=deltap)) +
    # Show all points
    # Show all points
    geom_point(aes(color=as.factor(Chr)), alpha=0.5, size=0.2) +
    scale_color_manual(values = rep(c("grey90", "grey75"), 22 )) +
    
    #geom_point(aes(colour=factor(seq(1, length(points.col)))), alpha=0.9, size=0.8) +
    #scale_color_manual(values = points.col) +
    # custom X axis:
    scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    # Add highlighted points
    geom_point(data=subset(plot_dat, is_highlight=="yes"), color=outlier_color, alpha=0.6, size=0.8) +
    facet_grid(cols = vars(Chr), scales = "free_x", space="free_x") +
    ylab(as.expression(bquote(.(x_label)~Delta~p)))+
    ylim(-0.5, 0.5)+
    theme_bw() +
    theme( 
      text=element_text(color ="black",size = 12),
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.border = element_blank(),
      axis.title.x=element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #panel.grid.major.y = element_blank(),
      #panel.grid.minor.y = element_blank()
    )
  #points.col <- colByValue(abs(don_sample$deltap), colors=c("white", "orange"))
}

manhattan_plot_bottom <- function(pname, outlier_name, outlier_color, x_label){
  #pname="ps_18_SGS_HC_NB.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat$SNP = paste0(dat$V1 , "_", dat$V2)
  colnames(dat) = c("Chr", "Pos", "AF", "BF", "deltap", "ps", "signif", "SNP")
  # A list of SNP of interest is provided
  #outlier_name="18_HC_NB_shared_outlier.txt"
  outlier_file <- read.delim(outlier_name, header = T, sep='\t')
  snpsOfInterest <- outlier_file$id
  
  don <- dat %>% 
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len=max(Pos)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(dat, ., by=c("Chr"="Chr")) %>%
    # Add a cumulative position of each SNP
    arrange(Chr, Pos) %>%
    mutate(BPcum=Pos+tot) %>%
    mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) 
  
  axisdf = don %>% group_by(Chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  sig_data <- don %>% 
    subset(is_highlight == "yes")
  notsig_data <- don %>% 
    subset(is_highlight == "no") %>%
    group_by(Chr) %>% 
    sample_frac(0.01)
  
  plot_dat <- bind_rows(sig_data, notsig_data)
  
  ggplot(plot_dat, aes(x=BPcum, y=deltap)) +
    # Show all points
    # Show all points
    geom_point(aes(color=as.factor(Chr)), alpha=0.5, size=0.2) +
    scale_color_manual(values = rep(c("grey90", "grey75"), 22 )) +
    
    #geom_point(aes(colour=factor(seq(1, length(points.col)))), alpha=0.9, size=0.8) +
    #scale_color_manual(values = points.col) +
    # custom X axis:
    scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    # Add highlighted points
    geom_point(data=subset(plot_dat, is_highlight=="yes"), color=outlier_color, alpha=0.6, size=0.8) +
    facet_grid(cols = vars(Chr), scales = "free_x", space="free_x") +
    ylab(as.expression(bquote(.(x_label)~Delta~p)))+
    xlab("Chromosome")+
    ylim(-0.5, 0.5)+
    theme_bw() +
    theme( 
      text=element_text(color ="black",size = 12),
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #panel.grid.major.y = element_blank(),
      #panel.grid.minor.y = element_blank()
    )
  #points.col <- colByValue(abs(don_sample$deltap), colors=c("white", "orange"))
}

mycolors <-  c("#231B12", "#9A9EAB", "#D0E1F9", "#DDBC95", "#B36687")
#  Sur_20     Sur_19      HC_21      HC_19     HC_18
#mycolor<-c("#A71B4B", "#E97302", "#4461A8", "#0BC0B3","#EAC728")
HC_NB_18 <- manhattan_plot("ps_18_SGS_HC_NB.txt", "18_HC_NB_shared_outlier.txt", "#231B12", "18 HC-NB")
HC_NB_19 <- manhattan_plot("ps_19_SGS_HC_NB.txt", "19_HC_NB_shared_outlier.txt", "#9A9EAB", "19 HC-NB")
HC_NB_21 <- manhattan_plot("ps_21_SGS_HC_NB.txt", "21_HC_NB_shared_outlier.txt", "#D0E1F9", "21 HC-NB")
Sur_Ref_19 <- manhattan_plot("ps_19_SGS_Sur_Ref.txt", "19_Sur_Ref_shared_outlier.txt", "#DDBC95", "19 Sur-Ref")
Sur_Ref_20 <- manhattan_plot_bottom("ps_20_SGS_Sur_Ref.txt", "20_Sur_Ref_shared_outlier.txt", "#B36687", "20 Sur-Ref")

#ylab <- HC_NB_18$labels$y
#HC_NB_18$labels$y <- HC_NB_19$labels$y <- HC_NB_21$labels$y <- Sur_Ref_19$labels$y <- Sur_Ref_20$labels$y <- " "

layout <-" 
AAAAAA
BBBBBB
CCCCCC
DDDDDD
EEEEEE
"
library(grid)
library(gtable)
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

p1 <- HC_NB_18 + HC_NB_19 + HC_NB_21 + Sur_Ref_19 + Sur_Ref_20 + plot_layout(design=layout, guides='collect') +
  theme(plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),units="inches"))
#grid::grid.draw(grid::textGrob("Chromosome", y=0.015, gp = gpar(col = "black", fontsize = 16)))
#grid::grid.draw(grid::textGrob(ylab, x=0.02, gp = gpar(col = "black", fontsize = 16)))

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_2")
tiff("figure2_part1.tiff", units="in", width=12, height=10, res=300)
p1
dev.off()
#graph2ppt(file="figure2_part1", width=10, height=12)
###################################
############ part 2 ###############
###################################
# part2 raincloud plot for local LD
library(ggpubr)
library(cowplot)
library(Rmisc)
library(export)
#install.packages("gghalves")
library(ggplot2)
library(gghalves)
library(patchwork)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/")

format <- function (pname1, pname2){
  #pname1 = "./15_LD_outlier_snps/Challenge_19.ngsld.output"
  dat1 = read.delim(pname1, header = FALSE, sep='\t')
  dat1 = dat1[,c(3,7)]
  dat1 = na.omit(dat1)
  dat1 = dat1[which(dat1$V3 < 500),]
  #pname2 = "./15_LD_random_snps/HC_19.97.ngsld.output"
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

comp_mean <- function(pop_name){
  random_mean <- c()
  for(i in seq(1,100)){
    df <- format(paste0("./15_LD_outlier_snps/", pop_name, ".ngsld.output"), paste0("./15_LD_random_snps/", pop_name ,"_", i, ".ngsld.output"))
    df2 <- summarySE(df, measurevar="V7", groupvars=c("group"))                            
    message(i)
    print(df2$V7[2])
    random_mean = c(random_mean, df2$V7[2])
  }
  
  df<-format(paste0("./15_LD_outlier_snps/", pop_name, ".ngsld.output"), paste0("./15_LD_random_snps/", pop_name ,"_1.ngsld.output"))
  observed_mean  <- summarySE(df, measurevar="V7", groupvars=c("group"))$V7[1]
  random_mean <- as.data.frame(random_mean)
  random_mean_point <- mean(random_mean$random_mean)
  message("observed_mean: ", observed_mean)
  message("random_mean_point: ", random_mean_point)
  
  p <- sum(random_mean > observed_mean)/dim(random_mean)[1]
  message("p-value: ", p)
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
  
  print(hist_p)
  graph2ppt(file=paste0(pop_name, "_r2_comp"), width=4, height=3) 
  return(random_mean)
  
}


Wild_18 <- comp_mean("Wild_18")
Wild_18_observed_mean <- 0.16892831248846
Wild_18_random_mean_point <- 0.169439022382044
Wild_18_p <- 0.57

Wild_19 <- comp_mean("Wild_19")
Wild_19_observed_mean <- 0.172259434747324
Wild_19_random_mean_point <- 0.172203491839576
Wild_19_p <- 0.05

Wild_21 <- comp_mean("Wild_21")
Wild_21_observed_mean <- 0.168348209629814
Wild_21_random_mean_point <- 0.167803875317812
Wild_21_p <- 0.45

Challenge_19 <- comp_mean("Challenge_19")
Challenge_19_observed_mean <- 0.180397355117793
Challenge_19_random_mean_point <- 0.175967811012246
Challenge_19_p <- 0.01
  
Challenge_20 <- comp_mean("Challenge_20")
Challenge_20_observed_mean <- 0.166551389685253
Challenge_20_random_mean_point <- 0.169392772462925
Challenge_20_p <- 0.87

NB_18 <- comp_mean("NB_18")
NB_18_observed_mean <- 0.178543120075959
NB_18_random_mean_point <- 0.182357994094058
NB_18_p <- 0.97
NB_19 <- comp_mean("NB_19")
NB_19_observed_mean <- 0.182948187546109
NB_19_random_mean_point <- 0.18412565453885
NB_19_p <- 0.7
NB_21 <- comp_mean("NB_21")
NB_21_observed_mean <- 0.172951037829669
NB_21_random_mean_point <- 0.177135175122511
NB_21_p <- 0.95
Ref_19 <- comp_mean("Ref_19")
Ref_19_observed_mean <- 0.185841918380583
Ref_19_random_mean_point <- 0.185181412348559
Ref_19_p <- 0.39
Ref_20 <- comp_mean("Ref_20")
Ref_20_observed_mean <- 0.164516829106079
Ref_20_random_mean_point <- 0.169563750184561
Ref_20_p <- 0.93

HC_18 <- comp_mean("HC_18")
HC_18_observed_mean <- 0.181902400634167
HC_18_random_mean_point <- 0.179888406784559
HC_18_p <- 0.15
HC_19 <- comp_mean("HC_19")
HC_19_observed_mean <- 0.18910645698106
HC_19_random_mean_point <- 0.181051261768928
HC_19_p <- 0
HC_21 <- comp_mean("HC_21")
HC_21_observed_mean <- 0.181695596813964
HC_21_random_mean_point <- 0.176651046238946
HC_21_p <- 0.03
Sur_19 <- comp_mean("Sur_19")
Sur_19_observed_mean <- 0.194435692839424
Sur_19_random_mean_point <- 0.189270750471226
Sur_19_p <- 0
Sur_20 <- comp_mean("Sur_20")
Sur_20_observed_mean <- 0.182546108292584
Sur_20_random_mean_point <- 0.177611800015779
Sur_20_p <- 0.01

df_total = rbind(Wild_18, Wild_19, Wild_21, Challenge_19, Challenge_20)
df_total$group = as.factor(rep(c("Wild_18", "Wild_19", "Wild_21", "Challenge_19", "Challenge_20"), c(dim(Wild_18)[1], dim(Wild_19)[1], dim(Wild_21)[1], dim(Challenge_19)[1], dim(Challenge_20)[1])))
colnames(df_total) = c("random_mean", "group")
df_total$group = factor(df_total$group, levels = rev(c("Wild_18", "Wild_19", "Wild_21", "Challenge_19", "Challenge_20")))

df_total = rbind(HC_18, HC_19, HC_21, Sur_19, Sur_20)
df_total$group = rep(c("HC_18", "HC_19", "HC_21", "Sur_19", "Sur_20"), c(dim(HC_18)[1], dim(HC_19)[1], dim(HC_21)[1], dim(Sur_19)[1], dim(Sur_20)[1]))
colnames(df_total) = c("random_mean", "group")
df_total$group = factor(df_total$group, levels = c("Sur_20" , "Sur_19", "HC_21", "HC_19", "HC_18"))

df_total_control = rbind(NB_18, NB_19, NB_21, Ref_19, Ref_20)
df_total_control$group = rep(c("NB_18", "NB_19", "NB_21", "Ref_19", "Ref_20"), c(dim(NB_18)[1], dim(NB_19)[1], dim(NB_21)[1], dim(Ref_19)[1], dim(Ref_20)[1]))
colnames(df_total_control) = c("random_mean", "group")
df_total_control$group = factor(df_total_control$group, levels = c("Ref_20", "Ref_19", "NB_21", "NB_19", "NB_18"))

########
# plot #
########
mycolor <-  rev(c("#231B12", "#9A9EAB", "#D0E1F9", "#DDBC95", "#B36687"))

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
p2_obs<-ggplot(df_total, aes(x=group, 
                         y=random_mean,
                         fill=group,
                         color=group))+
  scale_color_manual(values=(mycolor))+scale_fill_manual(values=(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                        side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  #geom_point(aes(x = as.numeric(group)-0.1, 
  #               y = random_mean,
  #               color = group), 
  #           position = position_jitter(width =0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  ylim(0.16, 0.2)+
  labs(y=expression(paste(Mean~linkage~disequilibrium~(r^2))), x = NULL)+
  geom_segment(aes(x = 5.4, y = HC_18_observed_mean, xend = 5.1, yend = HC_18_observed_mean, linetype="p>0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 4.4, y = HC_19_observed_mean, xend = 4.1, yend = HC_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 3.4, y = HC_21_observed_mean, xend = 3.1, yend = HC_21_observed_mean, linetype="p<0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 2.4, y = Sur_19_observed_mean, xend = 2.1, yend = Sur_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 1.4, y = Sur_20_observed_mean, xend = 1.1, yend = Sur_20_observed_mean, linetype="p<0.05"), colour="red", size=0.8)+
  scale_linetype_manual(
    expression(paste(Observed~mean~(r^2))), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p2_obs

p2_control<-ggplot(df_total_control, aes(x=group, 
                         y=random_mean,
                         fill=group,
                         color=group))+
  scale_color_manual(values=(mycolor))+scale_fill_manual(values=(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  #geom_point(aes(x = as.numeric(group)-0.1, 
  #               y = random_mean,
  #               color = group), 
  #           position = position_jitter(width =0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  ylim(0.16, 0.2)+
  labs(y=expression(paste(Mean~linkage~disequilibrium~(r^2))), x = NULL)+
  geom_segment(aes(x = 5.4, y = NB_18_observed_mean, xend = 5.1, yend = NB_18_observed_mean, linetype="p>0.05"), colour="red", size=0.8,)+
  geom_segment(aes(x = 4.4, y = NB_19_observed_mean, xend = 4.1, yend = NB_19_observed_mean, linetype="p>0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 3.4, y = NB_21_observed_mean, xend = 3.1, yend = NB_21_observed_mean, linetype="p>0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 2.4, y = Ref_19_observed_mean, xend = 2.1, yend = Ref_19_observed_mean, linetype="p>0.05"), colour="red", size=0.8)+
  geom_segment(aes(x = 1.4, y = Ref_20_observed_mean, xend = 1.1, yend = Ref_20_observed_mean, linetype="p>0.05"), colour="red", size=0.8)+
  scale_linetype_manual(
    expression(paste(Observed~mean~(r^2))), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p2_control

ylab <- p2_obs$labels$y
p2_obs$labels$y <- p2_control$labels$y <- " "

layout <-" 
AAABBB
AAABBB
AAABBB
AAABBB
"
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_2")
tiff("figure2_part4.jpg", units="in", width=16, height=8, res=300)
p2_obs + p2_control + plot_layout(design=layout, guides='collect')
#grid::grid.draw(grid::textGrob(expression(paste(Mean ~ linkage ~ disequilibrium ~ (r^2))), y=0.015, gp = gpar(col = "black", fontsize = 16)))
dev.off()

graph2ppt(file="figure2_part4", width=10, height=12)

pg <- ggplot_build(p2) 
pg$data[1][[1]]


#一半小提琴图的参数调整：
#position：位置调整，这里将其向右水平移动0.1；
#side：显示哪一侧， "I"代表左侧，"R"代表右侧，默认"I"；
#adjust：调整带宽，这里设为1.2使宽带略变平滑；
#trim：小提琴图尾部的数据修整，默认为"T",表示将尾部修整到数据范围；"F"表示不修剪尾部；
p_all<-ggplot(df_total, aes(x=group, 
                             y=random_mean,
                             fill=group,
                             color=group))+
  scale_color_manual(values=(mycolor))+scale_fill_manual(values=(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  #geom_point(aes(x = as.numeric(group)-0.1, 
  #               y = random_mean,
  #               color = group), 
  #           position = position_jitter(width =0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  ylim(0.155, 0.185)+
  labs(y=expression(paste(Mean~linkage~disequilibrium~(r^2))), x = NULL)+
  geom_segment(aes(x = 5.4, y = Wild_18_observed_mean, xend = 5.1, yend = Wild_18_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.02, "npc"), type="closed"))+
  geom_segment(aes(x = 4.4, y = Wild_19_observed_mean, xend = 4.1, yend = Wild_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.02, "npc"), type="closed"))+
  geom_segment(aes(x = 3.4, y = Wild_21_observed_mean, xend = 3.1, yend = Wild_21_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.02, "npc"), type="closed"))+
  geom_segment(aes(x = 2.4, y = Challenge_19_observed_mean, xend = 2.1, yend = Challenge_19_observed_mean, linetype="p<0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.02, "npc"), type="closed"))+
  geom_segment(aes(x = 1.4, y = Challenge_20_observed_mean, xend = 1.1, yend = Challenge_20_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.02, "npc"), type="closed"))+
  scale_linetype_manual(
    expression(paste(Observed~mean~(r^2))), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p_all

###################################
############ part 3 ###############
###################################
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
mycolor <-  rev(c("#231B12", "#9A9EAB", "#D0E1F9", "#DDBC95", "#B36687"))


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
  # geom_point(aes(x = as.numeric(group)-0.1, 
  #                y = random_mean,
  #                color = group), 
  #            position = position_jitter(width = 0.01),size =0.2, shape = 20) + 
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
  # geom_point(aes(x = as.numeric(group)-0.1, 
  #                y = random_mean,
  #                color = group), 
  #            position = position_jitter(width = 0.01),size =0.2, shape = 20) + 
  geom_boxplot(outlier.shape = NA, 
               width =0.1,
               alpha=0.7)+
  labs(y=expression(paste(Mean~B[0*', MAF'])), x = NULL)+
  geom_segment(aes(x = 5.4, y = NB_18_observed_mean, xend = 5.1, yend = NB_18_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 4.4, y = NB_19_observed_mean, xend = 4.1, yend = NB_19_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 3.4, y = NB_21_observed_mean, xend = 3.1, yend = NB_21_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 2.4, y = Ref_19_observed_mean, xend = 2.1, yend = Ref_19_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = 1.4, y = Ref_20_observed_mean, xend = 1.1, yend = Ref_20_observed_mean, linetype="p>0.05"), colour="red", size=0.8, arrow = arrow(length = unit(0.3, "cm")))+
  scale_linetype_manual(
    expression(paste(Observed~mean~B[0*', MAF'])), 
    values = c("p<0.05" = "solid", "p>0.05" = "dotted"), 
  ) + 
  coord_flip()+
  guides(fill="none")+
  labs(colour = "Population") +
  mytheme

p3_random

ylab <- p3_obs$labels$y
p3_obs$labels$y <- p3_random$labels$y <- " "

layout <-" 
AAABBB
AAABBB
AAABBB
AAABBB
"
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_2")
tiff("figure2_part5.jpg", units="in", width=16, height=8, res=300)
p3_obs + p3_random + plot_layout(design=layout, guides='collect')
#grid::grid.draw(grid::textGrob(expression(paste(Mean ~ linkage ~ disequilibrium ~ (r^2))), y=0.015, gp = gpar(col = "black", fontsize = 16)))
dev.off()

graph2ppt(file="figure2_part5", width=10, height=12)

###################################
############ part 4 ###############
###################################
# part4 delta-p and starting allele frequency distribution
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_2/delta_p")

df_plot <-function(input) {
  df <- read.delim(input, header = TRUE, sep='\t')
  
  toPlot <-
    df %>%
    mutate(direction = ifelse(deltap > 0, "UP", "DOWN")) %>%
    mutate(abs_deltap = abs(deltap)) %>%
    mutate(bin = cut(abs_deltap, breaks=seq(from=0, to=0.5, by=0.02))) %$%
    
    #cut(density(bs[,i], bw = 0.05)$x, breaks=seq(from=-0.3, to=1.3, by=0.01))
    
    table(bin, direction) %>%
    as.data.frame() %>%
    mutate(plotVal = ifelse(direction == "DOWN"
                            , -1*Freq
                            , Freq))
  toPlot$group <- gsub('_shared_outlier.txt','',input)
  return(toPlot)
}

d1 <- df_plot("18_HC_NB_shared_outlier.txt")
d2 <- df_plot("19_HC_NB_shared_outlier.txt")
d3 <- df_plot("21_HC_NB_shared_outlier.txt")
d4 <- df_plot("19_Sur_Ref_shared_outlier.txt")
d5 <- df_plot("20_Sur_Ref_shared_outlier.txt")

df <- rbind(d1, d2, d3, d4, d5)

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

mytheme<-theme_bw()+
  theme(
    axis.text.x = element_text(colour = "black", size=10),
    axis.text.y = element_text(colour = "black", size=10),
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


mycolors <-  c("#231B12", "#CDCDC0", "#D0E1F9", "#DDBC95", "#B36687")

order = c("18_HC_NB", "19_HC_NB", "21_HC_NB", "19_Sur_Ref", "20_Sur_Ref")
df$group <-factor(df$group, levels=order)

df %>%
  ggplot(aes(x = bin, y = plotVal, fill = factor(group))) +
  geom_bar(position="dodge", stat="identity")+
  #geom_col()+
  scale_fill_manual(values=mycolors)+
  theme_classic()+
  ylim(-120,300)+
  #facet_grid(group~., scales="free_y")+
  geom_hline(aes(yintercept=0), linetype=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4))+
  xlab(expression("Allele frequency changes "*Delta~italic(p)))+
  ylab("Counts")

df %>%
  ggplot(aes(x = bin, y = plotVal, fill = factor(group))) +
  #geom_bar(position="dodge", stat="identity")+
  geom_col()+
  scale_fill_manual(values=mycolors)+
  theme_classic()+
  ylim(-120,300)+
  facet_grid(group~., scales="free_y")+
  geom_hline(aes(yintercept=0), linetype=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4))+
  xlab(expression("Allele frequency changes "*Delta~italic(p)))+
  ylab("Counts")

#############################
# starting allele frequency #
#############################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/Figures/Figure_2/delta_p")
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(grid)
library(gtable)
library(export)
library(dplyr)
#pname = "./ps_NB.txt"
#dat1 = read.delim(pname, header = TRUE, sep='\t')
########################################
# starting and ending allele frequency #
########################################

# old figure for allele frequency patterns
# df_format <-function(input){
#   #input = "18_HC_NB_shared_outlier.txt"
#   dat = read.delim(input, header = TRUE, sep='\t')
#   df1 = data.frame(dat$p2, "starting allele frequency")
#   colnames(df1) = c("ps", "class")
#   df2 = data.frame(dat$p1, "ending allele frequency")
#   colnames(df2) = c("ps", "class")
#   df_plot = rbind(df1,df2)
#   df_plot$group <- gsub('_shared_outlier.txt','',input)
#   return(df_plot)
# }
# 
# 
# d1 <- df_format("18_HC_NB_shared_outlier.txt")
# d2 <- df_format("19_HC_NB_shared_outlier.txt")
# d3 <- df_format("21_HC_NB_shared_outlier.txt")
# d4 <- df_format("19_Sur_Ref_shared_outlier.txt")
# d5 <- df_format("20_Sur_Ref_shared_outlier.txt")
# 
# d2$ps
# 
# df <- rbind(d1, d2, d3, d4, d5)
# order = c("18_HC_NB", "19_HC_NB", "21_HC_NB", "19_Sur_Ref", "20_Sur_Ref")
# df$group <-factor(df$group, levels=order)
# 
# head(df)
# ggplot(df, aes(x=ps, color=class)) +  
#   geom_density(alpha = 0.2, size = 1.2) +
#   #geom_vline(aes(xintercept=grp.mean, color=class), linetype="dashed")+
#   scale_color_manual(values=c("#D09683", "#2D4262")) +
#   facet_grid(group~., scales="free_y")+
#   ylim(0,6)+
#   ylab("Density") + xlab("Allele frequency") +
#   theme_classic()+
#   theme(legend.position = "top")
########################################
# starting, ending, and neutral plots ##
########################################

# df_format <-function(input){
#   #input = "18_HC_NB_shared_outlier.txt"
#   dat = read.delim(input, header = TRUE, sep='\t')
#   df1 = data.frame(dat$p2, "starting allele frequency")
#   colnames(df1) = c("p2", "class")
#   df2 = data.frame(dat$p1, "ending allele frequency")
#   colnames(df2) = c("p1", "class")
#   df_plot = data.frame(df1,df2)
#   df_plot$group <- gsub('_shared_outlier.txt','',input)
#   return(df_plot)
# }
# 
# 
# d1 <- df_format("18_HC_NB_shared_outlier.txt")
# d2 <- df_format("19_HC_NB_shared_outlier.txt")
# d3 <- df_format("21_HC_NB_shared_outlier.txt")
# d4 <- df_format("19_Sur_Ref_shared_outlier.txt")
# d5 <- df_format("20_Sur_Ref_shared_outlier.txt")

# cor.test(d1$p2, d1$p1,
#          method = "pearson",
#          conf.level = 0.95)
# 
# # starting allele frequency
# plot(density(d1$p2, bw = 0.05), lwd = 2,
#      col = "red", main = "starting af")
# 
# # ending allele frequency
# plot(density(d1$p1, bw = 0.05), lwd = 2,
#      col = "red", main = "ending af")



plot_af <- function(outlier_input, full_input){
  #outlier_input <- "18_HC_NB_shared_outlier.txt"
  outlier_df = read.delim(outlier_input, header = TRUE, sep='\t')
  colnames(outlier_df) = c("chr", "pos", "p1", "p0", "deltap", "p_value", "fdr", "id")
  #full_input <-  "ps_18_SGS_HC_NB.txt"
  dat_full <-  read.delim(full_input, header = FALSE, sep='\t')
  
  nrep <- 100
  nsamp <- dim(outlier_df)[1]
  
  bs_p0 <- matrix(nrow=nsamp, ncol=nrep)
  avg_rep_p0 <- c()
  for(i in 1:nrep){
    bs_p0[,i] <- sample(dat_full$V4, size=nsamp, replace=FALSE)
    avg_rep_p0 <- c(avg_rep_p0, bs_p0[,i])
  }
  out_p0 <- list()
  for (i in 1: ncol(bs_p0)){
    out_p0[[i]] <- data.frame(x=density(bs_p0[,i], bw = 0.05)$x,
                           y= density(bs_p0[,i], bw = 0.05)$y,
                           bin = cut(density(bs_p0[,i], bw = 0.05)$x, breaks=seq(from=-0.3, to=1.3, by=0.01)))
  }
  out.new_p0 <- do.call(rbind, out_p0)
  densities.qtiles_p0 <- out.new_p0 %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(q05 = quantile(y, 0.025),
                     q50 = quantile(y, 0.5),
                     q95 = quantile(y, 0.975))
  densities.qtiles_p0$x <- seq(from=-0.1, to=-0.1+(dim(densities.qtiles_p0)[1]-1)*0.01, by=0.01)
  avg_perm_p0 <- c()
  for (i in 1: ncol(bs_p0)){
    avg_perm_p0 <- c(avg_perm_p0, bs_p0[,i])
  }
  
  bs_p1 <- matrix(nrow=nsamp, ncol=nrep)
  avg_rep_p1 <- c()
  for(i in 1:nrep){
    bs_p1[,i] <- sample(dat_full$V3, size=nsamp, replace=FALSE)
    avg_rep_p1 <- c(avg_rep_p1, bs_p1[,i])
  }
  out_p1 <- list()
  for (i in 1: ncol(bs_p1)){
    
    out_p1[[i]] <- data.frame(x=density(bs_p1[,i], bw = 0.05)$x,
                              y= density(bs_p1[,i], bw = 0.05)$y,
                              bin = cut(density(bs_p1[,i], bw = 0.05)$x, breaks=seq(from=-0.3, to=1.3, by=0.01)))
  }
  
  out.new_p1 <- do.call(rbind, out_p1)
  
  densities.qtiles_p1 <- out.new_p1 %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(q05 = quantile(y, 0.025),
                     q50 = quantile(y, 0.5),
                     q95 = quantile(y, 0.975))
  
  densities.qtiles_p1$x <- seq(from=-0.1, to=-0.1+(dim(densities.qtiles_p1)[1]-1)*0.01, by=0.01)
  
  avg_perm_p1 <- c()
  for (i in 1: ncol(bs_p1)){
    avg_perm_p1 <- c(avg_perm_p1, bs_p1[,i])
  }
  
  #png("~/urchin_af/figures/Fig_S04_unfolded_af.png", height=85, width=85, units="mm", res=300)
  plot(density(0:1), ylim=c(0,4.1),xlim=c(0,1), lwd=0,
       main="",
       ylab="",
       xlab="",
       cex.lab=1.1, cex.axis=1,
       xaxt="n",yaxt="n")
  
  axis(1, mgp=c(1.8, .2, 0), cex.axis=1,tcl=-0.2) # second is tick mark labels
  axis(2, mgp=c(1.8, .4, 0), cex.axis=1, tcl=-0.2)
  title(xlab="Allele frequency", line=1.5, cex.lab=1.2)
  title(ylab="Density", line=1.5, cex.lab=1.2)
  axis(2, mgp=c(1.8, .4, 0), cex.axis=1, tcl=-0.2)
  
  polygon(x=c(densities.qtiles_p0$x,rev(densities.qtiles_p0$x)),
          y=c(densities.qtiles_p0$q05,rev(densities.qtiles_p0$q95)),
          col=alpha("#5BC8AC", alpha=0.4),border=NA)
  lines(densities.qtiles_p0$x,densities.qtiles_p0$q50, col="#5BC8AC", lwd=3 )
  polygon(x=c(densities.qtiles_p1$x,rev(densities.qtiles_p1$x)),
          y=c(densities.qtiles_p1$q05,rev(densities.qtiles_p1$q95)),
          col=alpha("#F18D91", alpha=0.4),border=NA)
  lines(densities.qtiles_p1$x,densities.qtiles_p1$q50, col="#F18D91", lwd=3 )

  #lines(density(mydata$af_out), col=alpha("black", 1), lwd=3)
  lines(density(outlier_df$p1, bw=0.05), col=alpha("#CB0000", 1), lwd=3)
  lines(density(outlier_df$p0, bw=0.05), col=alpha("#1E656D", 1), lwd=3)
  #lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)
  abline(v=mean(avg_perm_p0), lty=2, col= "#5BC8AC", lwd=3) # mean: 0.459
  abline(v=mean(avg_perm_p1), lty=2, col= "#F18D91", lwd=3) # mean: 0.459
  abline(v=mean(outlier_df$p1), lty=2, col= "#CB0000", lwd=3) # mean: 0.248  note, shifting this slightly so visible in plot.
  abline(v=mean(outlier_df$p0), lty=2, col= "#1E656D", lwd=3) # mean: 0.251
  #abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3) # mean: 0.32
}
# mar(c(bottom, left, top, and right)) The default is c(5.1, 4.1, 4.1, 2.1).

#png("./Fig2_part4.png", height=1000, width=500, units="mm", res=300)
par(mfrow = c(5, 1), mar=c(1.2, 4.1, 4.1, 2.1))
plot_af("18_HC_NB_shared_outlier.txt", "ps_18_SGS_HC_NB.txt")
plot_af("19_HC_NB_shared_outlier.txt", "ps_19_SGS_HC_NB.txt")
plot_af("21_HC_NB_shared_outlier.txt", "ps_21_SGS_HC_NB.txt")
plot_af("19_Sur_Ref_shared_outlier.txt", "ps_19_SGS_Sur_Ref.txt")
plot_af("20_Sur_Ref_shared_outlier.txt", "ps_20_SGS_Sur_Ref.txt")
graph2ppt(file="Figure_2_part4",width=10,height=18)
#dev.off()
legend("right", c("Low salinity", "Relatively high salinity", "Neutral in low salinity", "Neutral in relatively high salinity"),
       col=c("#CB0000", "#1E656D", "royalblue3", "darkorchid2"), lty=1,
       cex=0.59, lwd=1.9)



mtext(text=bquote(paste('(',italic('b'),')')),
      side=3, line=0,
      cex=1.5,
      at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

dev.off()

#################################################################################################

#########
#
# unfolded maf
#
#########


mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)
cut_off <- 0.05/(9828)

mydata$D1_8_maf <- (sapply(mydata$D1_8_af,function(x)
  ifelse(x > 0.5, (1-x), x)))

selected_7 <- mydata[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off),]
selected_8 <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off),]
selected_both <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off),]
neutral <- mydata[which(mydata$pH8_selection_pval > cut_off & mydata$pH7_selection_pval > cut_off),]

# pull out allele freqs

selected_7_D7_8 <- selected_7$D7_8_af
selected_7_D7_7 <- selected_7$D7_7_af
selected_7_D1_8 <- selected_7$D1_8_af

selected_8_D7_8 <- selected_8$D7_8_af
selected_8_D7_7 <- selected_8$D7_7_af
selected_8_D1_8 <- selected_8$D1_8_af

selected_both_D7_8 <- selected_both$D7_8_af
selected_both_D7_7 <- selected_both$D7_7_af
selected_both_D1_8 <- selected_both$D1_8_af

selected_all <- unique(rbind(selected_7,selected_8,selected_both))
selected_all_D7_8 <- selected_all$D7_8_af
selected_all_D7_7 <- selected_all$D7_7_af
selected_all_D1_8 <- selected_all$D1_8_af

d7_8 <- (selected_all_D7_8 - selected_all_D1_8)
d7_7 <- (selected_all_D7_7 - selected_all_D1_8)

d7_8_both <- (selected_both_D7_8 - selected_both_D1_8)
d7_7_both <- (selected_both_D7_7 - selected_both_D1_8)

d7_8_s7 <- (selected_7_D7_8 - selected_7_D1_8)
d7_7_s7 <- (selected_7_D7_7 - selected_7_D1_8)

d7_8_s8 <- (selected_8_D7_8 - selected_8_D1_8)
d7_7_s8 <- (selected_8_D7_7 - selected_8_D1_8)


cor.test(d7_7, d7_8,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s7, d7_8_s7,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s8, d7_8_s8,
         method = "pearson",
         conf.level = 0.95)

cor.test(d7_7_both, d7_8_both,
         method = "pearson",
         conf.level = 0.95)


# corr genome wide:

d7_8_all <- (neutral$D7_8_af - neutral$D1_8_af)
d7_7_all <- (neutral$D7_7_af - neutral$D1_8_af)

cor.test(d7_8_all, d7_7_all)


# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_maf[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_maf[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_maf[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

mean(snp.sel_75)
mean(snp.sel_80)
mean(snp.sel_both)
mean(mydata$D1_8_maf)

sd(snp.sel_75)/sqrt(length(snp.sel_75))
sd(snp.sel_80)/sqrt(length(snp.sel_80))
sd(snp.sel_both)/sqrt(length(snp.sel_both))
sd(mydata$D1_8_maf)/sqrt(length(mydata$D1_8_maf))

nrep <- 1000
nsamp <- length(snp.sel_75)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
  bs[,i] <- sample(mydata$D1_8_maf, size=nsamp, replace=FALSE)
  avg_rep <- c(avg_rep, bs[,i])
}

# calc ks tests

p_75_ks <- c()
p_80_ks <- c()
p_both_ks <- c()

for (i in 1:ncol(bs)){
  p_75_ks[i] <- ks.test(snp.sel_75, bs[,i])$p.value
  p_80_ks[i] <- ks.test(snp.sel_80, bs[,i])$p.value
  p_both_ks[i] <- ks.test(snp.sel_both, bs[,i])$p.value
}

ks.test(snp.sel_75, mydata$D1_8_maf)
ks.test(snp.sel_80, mydata$D1_8_maf)
ks.test(snp.sel_both, mydata$D1_8_maf)
ks.test(snp.sel_both, snp.sel_80)
ks.test(snp.sel_both, snp.sel_75)
ks.test(snp.sel_75, snp.sel_80)

length(which(p_75_ks < 0.05))
length(which(p_80_ks < 0.05))
length(which(p_both_ks < 0.05))





