library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)

#################################
######## Fst violin plot ########
#################################
data_process <- function(headname, pop){
  # headname = "CS_UMFS_noinvers."
  name = paste0(headname, "150.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  DT[,9][DT[,9]<0] = 0 #@chnage
  dat <- data.frame(fst=DT[,9], Contrast=pop) # @change
  dat$fst <- as.numeric(dat$fst)
  return(dat)
}
# wild-domestic
#dat1 <- data_process("CS_HC_noinvers.", "CS_HC")
dat2 <- data_process("MEW_MES_noinvers.", "MEW_MES")
dat3 <- data_process("LIW_LIS_noinvers.", "LIW_LIS")
dat4 <- data_process("DBW_DBS_noinvers.", "DBW_DBS")
dat5 <- data_process("NCW_NCS_noinvers.","NCW_NCS")

plotdat = rbind(dat2, dat3, dat4, dat5)
plotdat$Contrast = factor(plotdat$Contrast, levels=c('MEW_MES', 'LIW_LIS', 'DBW_DBS', 'NCW_NCS' ))
jpeg("Fst_violin_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Contrast, y=fst, fill=Contrast)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Contrast", y = expression(F[italic(ST)]*'(150 SNPs/non-overlapping window)')) 
p +  scale_fill_manual(values=c("#87b374", "#FFD73B", "#FF9636", "#FF5C4D","#BFD7ED","#60A3D9","#0074B7","#003B73", "#A3EBB1","#18A558","#21B6A8","#116530")) + 
      scale_y_continuous(limits=c(0, 0.25)) +
      theme(legend.title = element_text(color = "Black", size = 9),
            legend.text = element_text(color = "Black", size = 9))
dev.off()

library(export)
graph2ppt(file="Fst_violin",width=8,height=3)
