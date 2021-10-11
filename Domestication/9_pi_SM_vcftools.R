library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
setwd("~/Documents/HG/Domestication/09_pi_SM_vcftools")
################################################
######## Pi violin plot for SM approach ########
################################################
data_left <- function(headname, pop){
  #headname = "MEW_MES_noinvers."
  name = paste0(headname, "150.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(Pi=DT[,6], Population=pop) # @change
  dat$Pi <- as.numeric(dat$Pi)
  return(dat)
}

data_right <- function(headname, pop){
  name = paste0(headname, "150.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(Pi=DT[,7], Population=pop) # @change
  dat$Pi <- as.numeric(dat$Pi)
  return(dat)
}

dat1 <- data_left("MEW_MES_noinvers.", "MEW")
dat2 <- data_left("LIW_LIS_noinvers.", "LIW")
dat3 <- data_left("DBW_DBS_noinvers.", "DBW")
dat4 <- data_left("NCW_NCS_noinvers.", "NCW")
dat5 <- data_right("MEW_MES_noinvers.", "MES")
dat6 <- data_right("LIW_LIS_noinvers.", "LIS")
dat7 <- data_right("DBW_DBS_noinvers.", "DBS")
dat8 <- data_right("NCW_NCS_noinvers.", "NCS")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
plotdat$Population = factor(plotdat$Population, levels=c('MEW', 'LIW', 'DBW', 'NCW','MES', 'LIS', 'DBS', 'NCS'  ))
jpeg("Pi_150snp_per_window_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=Pi, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = expression(pi*' (150 SNPs/non-overlapping window)')) 
p +  scale_fill_manual(values=c("#003f5c", "#2f4b7c", "#665191", "#a05195",
                                "#d45087","#f95d6a","#ff7c43","#ffa600")) + 
  scale_y_continuous(limits=c(0, 0.4)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()

library(export)
graph2ppt(file="Fst_violin",width=8,height=3)

################################################
######## Pi violin plot for vcftools    ########
################################################

data_process <- function(headname){
  #headname = "DBS"
  name = paste0(headname, ".windowed.pi")
  DT = read.delim(name, header = TRUE, sep='\t')
  colnames(DT) <- c("chr","st", "ed", "nsites", "pi" )
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(pi=DT$pi, Population=headname) # @change
  dat$pi <- as.numeric(dat$pi)
  return(dat)
}
# wild-domestic
dat1 <- data_process("MEW")
dat2 <- data_process("LIW")
dat3 <- data_process("DBW")
dat4 <- data_process("NCW")
dat5 <- data_process("MES")
dat6 <- data_process("LIS")
dat7 <- data_process("DBS")
dat8 <- data_process("NCS")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
plotdat$Population = factor(plotdat$Population, levels=c('MEW', 'LIW', 'DBW', 'NCW','MES', 'LIS', 'DBS', 'NCS'  ))
jpeg("pi_vcftools_ind_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=pi, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = expression(pi*' (500K/non-overlapping window)')) 
p +  scale_fill_manual(values=c("#003f5c", "#2f4b7c", "#665191", "#a05195",
                                "#d45087","#f95d6a","#ff7c43","#ffa600")) + 
  scale_y_continuous(limits=c(0, 0.00025)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()
