#################################
######## phenotype plot   #######
#################################
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/26_slim")

headname = "post_migration_juvenile_phenotype.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, Pheno=DT$V4, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$Pheno)
head(dat)

#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv <- ggplot(dat, aes(x = Gen, y = Pheno)) + 
  geom_line(aes(color = Population)) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-migration juvenile phenotypes")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))+
  geom_hline(yintercept = -9, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -0, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 9, linetype="dotted",  color = "red", size=0.5)


headname = "post_selection_adult_phenotype.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, Pheno=DT$V4, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$Pheno)
head(dat)

#jpeg("post_selection_adult_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
adult <- ggplot(dat, aes(x = Gen, y = Pheno)) + 
  geom_line(aes(color = Population)) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-selection adult phenotypes")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))+
  geom_hline(yintercept = -9, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -0, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 9, linetype="dotted",  color = "red", size=0.5)

jpeg("phenotype_plot.jpg", width = 8, height = 5, units = 'in', res = 150)
ggarrange(juv, adult, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)
dev.off()

########## boxplot ###########

#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv <- ggplot(dat, aes(x = Gen, y = Pheno)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-migration juvenile phenotypes")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))+
  geom_hline(yintercept = -9, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -0, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 9, linetype="dotted",  color = "red", size=0.5)


headname = "post_selection_adult_phenotype.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, Pheno=DT$V4, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$Pheno)
head(dat)

#jpeg("post_selection_adult_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
adult <- ggplot(dat, aes(x = Gen, y = Pheno)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-selection adult phenotypes")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))+
  geom_hline(yintercept = -9, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = -0, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 3, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 6, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 9, linetype="dotted",  color = "red", size=0.5)

jpeg("phenotype_box_plot.jpg", width = 8, height = 5, units = 'in', res = 150)
ggarrange(juv, adult, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)
dev.off()
