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
    select(-chr_len) %>%
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
    scale_color_manual(values = rep(c("grey90", "grey50"), 22 )) +
    
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
    select(-chr_len) %>%
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
    scale_color_manual(values = rep(c("grey90", "grey50"), 22 )) +
    
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

         #  Sur_20     Sur_19      HC_21      HC_19     HC_18
#mycolor<-c("#A71B4B", "#E97302", "#4461A8", "#0BC0B3","#EAC728")
HC_NB_18 <- manhattan_plot("ps_18_SGS_HC_NB.txt", "18_HC_NB_shared_outlier.txt", "#EAC728", "18 HC-NB")
HC_NB_19 <- manhattan_plot("ps_19_SGS_HC_NB.txt", "19_HC_NB_shared_outlier.txt", "#0BC0B3", "19 HC-NB")
HC_NB_21 <- manhattan_plot("ps_21_SGS_HC_NB.txt", "21_HC_NB_shared_outlier.txt", "#4461A8", "21 HC-NB")
Sur_Ref_19 <- manhattan_plot("ps_19_SGS_Sur_Ref.txt", "19_Sur_Ref_shared_outlier.txt", "#E97302", "19 Sur-Ref")
Sur_Ref_20 <- manhattan_plot_bottom("ps_20_SGS_Sur_Ref.txt", "20_Sur_Ref_shared_outlier.txt", "#A71B4B", "20 Sur-Ref")

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
p1
