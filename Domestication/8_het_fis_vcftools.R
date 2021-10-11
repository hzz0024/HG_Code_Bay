library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)

#################################
######## Het violin plot ########
#################################
data_process <- function(headname){
  #headname = "DBS"
  name = paste0(headname, ".het")
  DT = read.delim(name, header = TRUE, sep='\t')
  colnames(DT) <- c("ind","ho", "he", "nsites", "f" )
  DT$het = DT$ho/DT$nsites
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(het=DT$het, Population=headname) # @change
  dat$het <- as.numeric(dat$het)
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
jpeg("Het_ind_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=het, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = 'Heterozygosity') 
p +  scale_fill_manual(values=c("#003f5c", "#2f4b7c", "#665191", "#a05195",
                                "#d45087","#f95d6a","#ff7c43","#ffa600")) + 
  scale_y_continuous(limits=c(0.65, 0.85)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()

library(export)
graph2ppt(file="Het_ind_plot",width=8,height=3)

#################################
######## Fis violin plot ########
#################################
data_process <- function(headname){
  #headname = "DBS"
  name = paste0(headname, ".het")
  DT = read.delim(name, header = TRUE, sep='\t')
  colnames(DT) <- c("ind","ho", "he", "nsites", "f" )
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(Fis=DT$f, Population=headname) # @change
  dat$Fis <- as.numeric(dat$Fis)
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
jpeg("Fis_ind_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=Fis, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = 'Inbreeding coefficient') 
p +  scale_fill_manual(values=c("#003f5c", "#2f4b7c", "#665191", "#a05195",
                                "#d45087","#f95d6a","#ff7c43","#ffa600")) + 
  scale_y_continuous(limits=c(0, 0.5)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()

library(export)
graph2ppt(file="Het_ind_plot",width=8,height=3)
