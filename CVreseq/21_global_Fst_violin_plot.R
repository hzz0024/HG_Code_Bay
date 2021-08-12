library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)

data_process <- function(headname, pop){
  # headname = "CS_UMFS_noinvers."
  name = paste0(headname, "150.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  DT[,9][DT[,9]<0] = 0.000001 #@chnage
  dat <- data.frame(fst=DT[,9], Contrast=pop) # @change
  dat$fst <- as.numeric(dat$fst)
  return(dat)
}

dat1 <- data_process("CS_HC_noinvers.", "CS_HC")
dat2 <- data_process("CLP_HCVA_noinvers.", "HCVA_CLP")
dat3 <- data_process("CL_SL_noinvers.", "CL_SL")
dat4 <- data_process("CS_UMFS_noinvers.", "CS_UMFS")
dat5 <- data_process("CS_DEBY_noinvers.","CS_DEBY")
dat6 <- data_process("CS_NEH_noinvers.", "CS_NEH")
dat7 <- data_process("NEH_UMFS_noinvers.", "NEH_UMFS")
dat8 <- data_process("DEBY_NEH_noinvers.", "NEH_DEBY")
dat9 <- data_process("CL_CS_noinvers.", "CS_CL")
dat10 <- data_process("NEH_OBOYS2_noinvers.", "NEH_OBOYS")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10)
mean1 = mean(rbind(dat1, dat2, dat3)[,1])
mean2 = mean(rbind(dat4, dat5, dat6)[,1])
mean3 = mean(rbind(dat7, dat8)[,1])
mean4 = mean(rbind(dat9, dat10)[,1])
plotdat$Contrast = factor(plotdat$Contrast, levels=c('CS_HC', 'CL_SL', 'HCVA_CLP', 'CS_DEBY', 'CS_NEH', 'CS_UMFS',"NEH_UMFS","NEH_DEBY", "CS_CL", "NEH_OBOYS"))
# 
#                      Population=c(rep('CS_HC',length(dat$fst)), 
#                                   rep('HCVA_CLP',length(dat$fst)),
#                                   rep('CL_SL',length(dat$fst)),
#                                   rep('CS_UMFS',length(dat$fst)),
#                                   rep('CS_DEBY',length(dat$fst)),
#                                   rep('CS_NEH',length(dat$fst))))

jpeg("Fst_box.jpg", width = 16, height = 9, units = 'in', res = 300)
p <- ggplot(plotdat, aes(x=Contrast, y=fst, fill=Contrast)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(width = 0.3, outlier.size = 0.05) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Contrast", y = "Fst (150 SNPs/non-overlapping window)") 
p +  scale_fill_manual(values=c("#FFD68A", "#FFB52E", "#FFA500", "#89CFF0","#38AEE6","#1B99D4","#5DE31D","#3B9212", "#FFC5D0", "#FB6090")) + 
     scale_y_continuous(limits=c(0, 1)) + 
     geom_hline(yintercept=mean1, linetype="solid",color = "#FFB52E")+ 
     geom_hline(yintercept=mean2, linetype="solid",color = "#38AEE6")+ 
     geom_hline(yintercept=mean3, linetype="solid",color = "#5DE31D")+ 
     geom_hline(yintercept=mean4, linetype="solid",color = "#FFC5D0")
graph2ppt(file="Fst_box",width=12,height=4)
dev.off()

############################ final run ############################ 
data_process <- function(headname, pop){
  # headname = "CS_UMFS_noinvers."
  name = paste0(headname, "150.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  DT[,9][DT[,9]<0] = 0.000001 #@chnage
  dat <- data.frame(fst=DT[,9], Contrast=pop) # @change
  dat$fst <- as.numeric(dat$fst)
  return(dat)
}
# wild-wild
dat1 <- data_process("CS_HC_noinvers.", "CS_HC")
dat2 <- data_process("CL_SL_noinvers.", "CL_SL")
dat3 <- data_process("HI_SM_noinvers.", "HI_SM")
dat4 <- data_process("CL_CS_noinvers.", "CL_CS")
# wild-domestic
dat5 <- data_process("CS_DEBY_noinvers.","CS_DEBY")
dat6 <- data_process("CL_OBOYS2_noinvers.", "CL_OBOYS2")
dat7 <- data_process("HI_UMFS_noinvers.", "CS_UMFS")
dat8 <- data_process("CLP_LOLA_noinvers.", "CS_NEH")
# domestic-domestic
dat9 <- data_process("LOLA_OBOYS2_noinvers.", "LOLA_OBOYS2")
dat10 <- data_process("NEH_UMFS_noinvers.", "NEH_UMFS")
dat11 <- data_process("DEBY_NEH_noinvers.", "NEH_DEBY")
dat12 <- data_process("OBOYS2_UMFS_noinvers.", "OBOYS2_UMFS")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12)
plotdat$Contrast = factor(plotdat$Contrast, levels=c('CS_HC', 'CL_SL', 'HI_SM', 'CL_CS', 'CS_DEBY', "CL_OBOYS2","CS_UMFS",'CS_NEH',"LOLA_OBOYS2","NEH_UMFS","NEH_DEBY","OBOYS2_UMFS" ))
# 
#                      Population=c(rep('CS_HC',length(dat$fst)), 
#                                   rep('HCVA_CLP',length(dat$fst)),
#                                   rep('CL_SL',length(dat$fst)),
#                                   rep('CS_UMFS',length(dat$fst)),
#                                   rep('CS_DEBY',length(dat$fst)),
#                                   rep('CS_NEH',length(dat$fst))))

jpeg("Fst_box.jpg", width = 5, height = 2.5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Contrast, y=fst, fill=Contrast)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(width = 0.3, outlier.size = 0.05) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Contrast", y = "Fst (150 SNPs/non-overlapping window)") 
p +  scale_fill_manual(values=c("#FFED01", "#FFD73B", "#FF9636", "#FF5C4D","#BFD7ED","#60A3D9","#0074B7","#003B73", "#A3EBB1","#18A558","#21B6A8","#116530")) + 
      scale_y_continuous(limits=c(0, 1)) +
      theme(legend.title = element_text(color = "Black", size = 9),
            legend.text = element_text(color = "Black", size = 9))
dev.off()

library(export)
graph2ppt(file="Fst_violin",width=8,height=3)
