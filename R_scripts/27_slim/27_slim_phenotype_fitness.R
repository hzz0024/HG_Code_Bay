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
dat = dat[dat$Population != "Pop1" & dat$Population != "Pop8", ]

#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv <- ggplot(dat, aes(x = Gen, y = Pheno)) + 
  geom_line(aes(color = Population)) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-migration juvenile phenotypes")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))+
  geom_hline(yintercept = -7, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 7, linetype="dotted",  color = "red", size=0.5)
  
headname = "post_selection_adult_phenotype_fitness.txt"
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
  geom_hline(yintercept = -7, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 7, linetype="dotted",  color = "red", size=0.5)


#jpeg("phenotype_plot.jpg", width = 8, height = 5, units = 'in', res = 150)
ggarrange(juv, adult, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)
#dev.off()

########## boxplot ###########

headname = "post_migration_juvenile_phenotype.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, Pheno=DT$V4, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$Pheno)
head(dat)
dat = dat[dat$Population != "Pop1" & dat$Population != "Pop8", ]
dat = dat[dat$Gen > 500, ]
head(dat)
#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv_phe <- ggplot(dat, aes(x = Population, y = Pheno)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x=NULL, y = "Post-migration juvenile phenotypes")+
  theme(legend.position = "None")+
  geom_hline(yintercept = -7, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 7, linetype="dotted",  color = "red", size=0.5)

headname = "post_selection_adult_phenotype_fitness.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, Pheno=DT$V4, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$Pheno)
dat = dat[dat$Gen > 500, ]
head(dat)

#jpeg("post_selection_adult_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
adult_phe <- ggplot(dat, aes(x = Population, y = Pheno)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x=NULL, y = "Post-selection adult phenotypes")+
  theme(legend.position = "None")+
  geom_hline(yintercept = -7, linetype="dotted",  color = "red", size=0.5)+
  geom_hline(yintercept = 7, linetype="dotted",  color = "red", size=0.5)

jpeg("phenotype_box_plot.jpg", width = 8, height = 3, units = 'in', res = 150)
ggarrange(juv_phe, adult_phe, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)
dev.off()

###############################################################################
###############################################################################
#rm(list=ls())
headname = "post_selection_adult_phenotype_fitness.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, fit=DT$V6, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$fit)
dat=dat[dat$Gen> 500,]
head(dat)

#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
adult <- ggplot(dat, aes(x = Gen, y = fit)) + 
  geom_line(aes(color = Population)) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-selection adult fitness")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))

headname = "post_selection_adult_phenotype_fitness.txt"
DT = read.delim(headname, header = FALSE, sep=',')
headname1 = "post_migration_juvenile_fitness.txt"
DT1 = read.delim(headname1, header = FALSE, sep=',')

all_vals = c()
for(row in seq(dim(DT1)[1])){
  vals = c(t(DT1[row,2:7]))
  all_vals = c(all_vals, vals)
}
DT$V7 = all_vals

dat <- data.frame(Gen=DT$V1, fit=DT$V7, Population=DT$V2) # @change
dat$fit <- as.numeric(dat$fit)
dat=dat[dat$Gen > 500,]
head(dat)

#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv <- ggplot(dat, aes(x = Gen, y = fit)) + 
  geom_line(aes(color = Population)) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x="Generations", y = "Post-migration juvenile fitness")+
  theme(legend.title = element_text(color = "Black", size = 9),
        legend.text = element_text(color = "Black", size = 9))

#jpeg("fitness_time_plot.jpg", width = 8, height = 5, units = 'in', res = 150)
ggarrange(juv, adult, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)
#dev.off()


######################## box plot ##########################

#rm(list=ls())
headname = "post_selection_adult_phenotype_fitness.txt"
DT = read.delim(headname, header = FALSE, sep=',')
dat <- data.frame(Gen=DT$V1, fit=DT$V6, Population=DT$V2) # @change
dat$Pheno <- as.numeric(dat$fit)
dat=dat[dat$Gen > 500,]
head(dat)

#jpeg("post_selection_adult_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
adult_fit <- ggplot(dat, aes(x = Population, y = fit)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x=NULL, y = "Post-selection adult fitness")+
  theme(legend.position = "None")

headname = "post_selection_adult_phenotype_fitness.txt"
DT = read.delim(headname, header = FALSE, sep=',')
headname1 = "post_migration_juvenile_fitness.txt"
DT1 = read.delim(headname1, header = FALSE, sep=',')

all_vals = c()
for(row in seq(dim(DT1)[1])){
  vals = c(t(DT1[row,2:7]))
  all_vals = c(all_vals, vals)
}
DT$V7 = all_vals

dat <- data.frame(Gen=DT$V1, fit=DT$V7, Population=DT$V2) # @change
dat$fit <- as.numeric(dat$fit)
dat=dat[dat$Gen > 500, ]
head(dat)


#jpeg("post_migration_juvenile_phenotype_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
juv_fit <- ggplot(dat, aes(x = Population, y = fit)) + 
  geom_boxplot(aes(color = Population), outlier.alpha = 0.001) +
  scale_color_manual(values = c("#10588f", "#737aac","#b1a1c7", "#e5cde4","#e3a0c4", "#e16f91", "#de425b")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x=NULL, y = "Post-migration juvenile fitness")+
  theme(legend.position = "None")


#jpeg("fitness_box_plot.jpg", width = 8, height = 5, units = 'in', res = 150)
ggarrange(juv_fit, adult_fit, 
          labels = c("A", "B"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)

jpeg("pheno_fitness_box_plot.jpg", width = 14, height = 8, units = 'in', res = 300)
ggarrange(juv_phe, adult_phe, juv_fit, adult_fit, 
          labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 2)
dev.off()
