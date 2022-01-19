library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/07_heterozygosity_fis_vcftools")
#################################
######## Het violin plot ########
#################################
data_process <- function(headname){
  #headname = "DBS"
  name = paste0(headname, ".het")
  DT = read.delim(name, header = TRUE, sep='\t')
  colnames(DT) <- c("ind","ho", "hoe", "nsites", "f" )
  DT$het = (DT$nsites - DT$ho)/DT$nsites
  DT <- DT[complete.cases(DT), ]
  dat <- data.frame(het=DT$het, Population=headname) # @change
  dat$het <- as.numeric(dat$het)
  return(dat)
}
# wild-domestic
dat1 <- data_process("MEW1")
dat2 <- data_process("MEW2")
dat3 <- data_process("LIW1")
dat4 <- data_process("LIW2")
dat5 <- data_process("DBW1")
dat6 <- data_process("DBW2")
dat7 <- data_process("NCW1")
dat8 <- data_process("NCW2")
dat9 <- data_process("UMFS")
dat10 <- data_process("MEH2")
dat11 <- data_process("NEH1")
dat12 <- data_process("NEH2")
dat13 <- data_process("DBX1")
dat14 <- data_process("DBX2")
dat15 <- data_process("DBX3")
dat16 <- data_process("UNC1")
dat17 <- data_process("UNC2")

#########   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         UMFS      MEH2       NEH1      NEH2        DBX1      DBX2         DBX3      UNC1       UNC2
col <- c( "#1BA3C6", "#2CB5C0",  "#30BCAD", "#21B087", "#33A65C", "#57A337", "#A2B627", "#F8B620", "#F8B620", "#F89217","#F06719", "#E03426",  "#F64971", "#FC719E", "#EB73B3", "#CE69BE", "#A26DC2")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12, dat13, dat14, dat15, dat16 ,dat17)
plotdat$Population = factor(plotdat$Population, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "UMFS", "MEH2", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
jpeg("Het_ind_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=het, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = 'Heterozygosity') 
p +  scale_fill_manual(values=col) + 
  scale_y_continuous(limits=c(0.2, 0.4)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()

p +  scale_fill_manual(values=col) + 
  scale_y_continuous(limits=c(0.2, 0.4)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
graph2ppt(file="Het_ind_plot",width=8,height=5)

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
dat1 <- data_process("MEW1")
dat2 <- data_process("MEW2")
dat3 <- data_process("LIW1")
dat4 <- data_process("LIW2")
dat5 <- data_process("DBW1")
dat6 <- data_process("DBW2")
dat7 <- data_process("NCW1")
dat8 <- data_process("NCW2")
dat9 <- data_process("UMFS")
dat10 <- data_process("MEH2")
dat11 <- data_process("NEH1")
dat12 <- data_process("NEH2")
dat13 <- data_process("DBX1")
dat14 <- data_process("DBX2")
dat15 <- data_process("DBX3")
dat16 <- data_process("UNC1")
dat17 <- data_process("UNC2")

plotdat = rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12, dat13, dat14, dat15, dat16 ,dat17)
plotdat$Population = factor(plotdat$Population, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "UMFS", "MEH2", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
jpeg("Fis_ind_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=Fis, fill=Population)) +
  #geom_violin(trim=FALSE, width = 1) + 
  geom_boxplot(notch=FALSE, width = 0.8, outlier.size = 0.1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = median, geom='point', shape=20, size=3, color='red', fill='red') +
  labs(x="Population", y = 'Inbreeding coefficient') 
p +  scale_fill_manual(values=col) + 
  scale_y_continuous(limits=c(-0.2, 0.2)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
dev.off()

p +  scale_fill_manual(values=col) + 
  scale_y_continuous(limits=c(-0.2, 0.2)) +
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))
graph2ppt(file="Fis_ind_plot",width=8,height=5)
