####################
# Examine the data #
####################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/01_pca")
# check the sample id
# note the sample_id_all_842.txt is created by extracting the samples from vcf file. see bcftools code below
# bcftools query -l genetyped_data_all_samples.vcf > sample_id_all_842.txt
# genetyped_data_all_samples.vcf is the direct output from Axiom Analysis Suite
file1 = "sample_id_all_842.txt"
dat1 = read.delim(file1, header = FALSE, sep='\t')
# count how many sample in the vcf 
paste0("number of samples in the vcf file is ", length(dat1$V1))
# check the sample id
head(dat1$V1)

# extract population information
# Here I used the "960_samples_QC_report.txt" prodvided by Zhenwei as a reference database. However, it only has 958 
# sample records (no individual DBX2-21 and MEW1-1). May need ask Zhenwei why it happens. 
require(dplyr)
file2 = "960_samples_QC_report.txt"
dat2 = read.delim(file2, header = TRUE, sep='\t')
head(dat2)
# it has a lot of columns but what we need is the matching sample id and population information
length(dat2$Sample.Filename[which(dat2$Sample.Filename %in% dat1$V1)])
# only extract the Sample.Filename, Rename, Ethnicity columns
array_dat <- dat2[which(dat2$Sample.Filename %in% dat1$V1),c(1,74,82)]
head(array_dat)
# count how many populations in the file - 28 here
apply(array_dat, 2, function(x) length(unique(x)))
# count how many individuals in each population
array_dat %>% count(Ethnicity)    

####################
#     PCA plots    #
####################


# for all sample n=842 
# vcf
vcf.fn <- "genetyped_data_all_maf05_maxmiss07.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_all_maf05_maxmiss07.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_all_maf05_maxmiss07.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_all_maf05_maxmiss07.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")

# plot(tab_pop$EV1,tab_pop$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""), 
#      pch=as.numeric(factor(population_nonum)), col=colors_chosen[as.numeric(factor(tab_pop$pop))])
# tmp = stringr::str_remove(levels(factor(tab_pop$pop)), "[0-9]+")
# levels(factor(tmp)) == levels(factor(population_nonum))
# legend("bottomright", legend=levels(factor(tab_pop$pop)), pch=as.numeric(factor(tmp)), col=colors_chosen, cex = 0.6, ncol=2)

# discrete color
library(RColorBrewer)
n <- 28
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
col_gradient = viridis_pal(option = "C")(28)  # n = number of colors seeked
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "CBW1", "CBW2", "NCW1", "NCW2","UMFS", "MEH2", "NYH1", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "VF16", "VF17", "VF18", "VF19", "VF20", "UNC1", "UNC2", "LGF")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop1, shape=pop1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = col_gradient)+
  scale_shape_manual(values=c(rep(15:18, 7)), breaks=order1)+
  #scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p1

order2 = c("MEW", "LIW", "DBW", "CBW", "NCW", "UMFS", "MEH", "NYH", "NEH", "DBX", "VF", "UNC", "LGF")
tab_pop$population_nonum1 <-factor(tab_pop$population_nonum, levels=order2)
# PC1-2 for locations
col_gradient = viridis_pal(option = "C")(13)
p2<-ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color = population_nonum1, shape=population_nonum1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values = c(col_gradient))+
  scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p2
grid.arrange(p1, p2, ncol=2)
require(gridExtra)
jpeg("PCA_all_842.jpg", width = 10, height = 8, units = 'in', res = 300)
p1
dev.off()

# for subsets n=606
# vcf
vcf.fn <- "genetyped_data_n_606_maf01_maxmiss07.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_606_maf01_maxmiss07.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_606_maf01_maxmiss07.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_606_maf01_maxmiss07.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_606.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_606_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

# plot(tab_pop$EV1,tab_pop$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""), 
#      pch=as.numeric(factor(population_nonum)), col=colors_chosen[as.numeric(factor(tab_pop$pop))])
# tmp = stringr::str_remove(levels(factor(tab_pop$pop)), "[0-9]+")
# levels(factor(tmp)) == levels(factor(population_nonum))
# legend("bottomright", legend=levels(factor(tab_pop$pop)), pch=as.numeric(factor(tmp)), col=colors_chosen, cex = 0.6, ncol=2)

# discrete color
library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
col_gradient = viridis_pal(option = "C")(20)  # n = number of colors seeked
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "CBW1", "CBW2", "NCW1", "NCW2","UMFS", "MEH2", "NYH1", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop1, shape=pop1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = col_gradient)+
  scale_shape_manual(values=c(rep(15:18, 5)), breaks=order1)+
  #scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p1

jpeg("PCA_606.jpg", width = 10, height = 8, units = 'in', res = 300)
p1
dev.off()

# for subsets n=576
# vcf
vcf.fn <- "genetyped_data_n_576_maf01_maxmiss07.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_576_maf01_maxmiss07.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_576_maf01_maxmiss07.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_576_maf01_maxmiss07.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_576.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_576_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

# plot(tab_pop$EV1,tab_pop$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""), 
#      pch=as.numeric(factor(population_nonum)), col=colors_chosen[as.numeric(factor(tab_pop$pop))])
# tmp = stringr::str_remove(levels(factor(tab_pop$pop)), "[0-9]+")
# levels(factor(tmp)) == levels(factor(population_nonum))
# legend("bottomright", legend=levels(factor(tab_pop$pop)), pch=as.numeric(factor(tmp)), col=colors_chosen, cex = 0.6, ncol=2)

# discrete color
library(RColorBrewer)
n <- 19
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
col_gradient = viridis_pal(option = "C")(19)  # n = number of colors seeked
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "UMFS", "MEH2", "LIW1", "LIW2", "NEH1", "NEH2", "DBW1", "DBW2", "DBX1", "DBX2", "DBX3", "CBW1", "CBW2", "NCW1", "NCW2", "UNC1", "UNC2")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop1, shape=pop1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = col_gradient)+
  scale_shape_manual(values=c(rep(15:18, 5)), breaks=order1)+
  #scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p1

jpeg("PCA_576.jpg", width = 10, height = 8, units = 'in', res = 300)
p1
dev.off()

# for subsets n=514
# vcf
vcf.fn <- "genetyped_data_n_514_maf01_maxmiss07.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_514_maf01_maxmiss07.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_514_maf01_maxmiss07.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_514_maf01_maxmiss07.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_514.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_514_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

# plot(tab_pop$EV1,tab_pop$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""), 
#      pch=as.numeric(factor(population_nonum)), col=colors_chosen[as.numeric(factor(tab_pop$pop))])
# tmp = stringr::str_remove(levels(factor(tab_pop$pop)), "[0-9]+")
# levels(factor(tmp)) == levels(factor(population_nonum))
# legend("bottomright", legend=levels(factor(tab_pop$pop)), pch=as.numeric(factor(tmp)), col=colors_chosen, cex = 0.6, ncol=2)

# discrete color
library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
library(dplyr) 
col_gradient = viridis_pal(option = "C")(17)  # n = number of colors seeked
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "UMFS", "MEH2", "LIW1", "LIW2", "NEH1", "NEH2", "DBW1", "DBW2", "DBX1", "DBX2", "DBX3", "NCW1", "NCW2", "UNC1", "UNC2")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop1, shape=pop1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = col_gradient)+
  scale_shape_manual(values=c(rep(15:18, 5)), breaks=order1)+
  #scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p1

jpeg("PCA_514.jpg", width = 10, height = 8, units = 'in', res = 300)
p1
dev.off()


# for subsets n=252, 220, 188
# vcf
vcf.fn <- "genetyped_data_66k.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_66k.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_66k.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_66k.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_66k.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_66k_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

# plot(tab_pop$EV1,tab_pop$EV2, cex=1.4,xlab=paste("set1 PC1 - ",round(pc.percent[1],3),"%",sep=""), ylab=paste("set1 PC2 - ",round(pc.percent[2],3),"%",sep=""), 
#      pch=as.numeric(factor(population_nonum)), col=colors_chosen[as.numeric(factor(tab_pop$pop))])
# tmp = stringr::str_remove(levels(factor(tab_pop$pop)), "[0-9]+")
# levels(factor(tmp)) == levels(factor(population_nonum))
# legend("bottomright", legend=levels(factor(tab_pop$pop)), pch=as.numeric(factor(tmp)), col=colors_chosen, cex = 0.6, ncol=2)

# discrete color
library(RColorBrewer)
n <- 8
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
#col_gradient = viridis_pal(option = "C")(9)  # n = number of colors seeked
#col_gradient <- c("#26497a", "#904994", "#ea426d", "#ff7d00", "#00bf0d", "#82de86", "#0080b8", "#86bde2")
col_gradient <- c("#26497a", "#904994", "#ea426d", "#ff7d00", "#00bf0d", "#82de86", "#86bde2")
#col_gradient <- c("#26497a", "#904994", "#ea426d", "#ff7d00", "#82de86", "#86bde2")
# PC1-2 for individual populations
order1 = c("DBW1", "DBW2", "LIW1", "LIW2", "NEH1", "NEH2", "DBX2")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop1, shape=pop1))+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = col_gradient)+
  scale_shape_manual(values=c(rep(15:18, 5)), breaks=order1)+
  #scale_shape_manual(values=c(rep(15:18, 7)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=15),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))+
  theme(text = element_text(size = 20)) 
  
p1

jpeg("PCA_66k.jpg", width = 10, height = 8, units = 'in', res = 300)
p1
dev.off()

#######################
##### formal run ######
#######################
# discrete color
library(RColorBrewer)
# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
library(dplyr) 
library(tibble)
library(patchwork)
library(export)
# see https://github.com/clairemerot/Intro2PopGenomics/blob/master/3.2.3/PCA/script_PCA_from_vcf.R for original code
library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

# this is the function for ellipse
geom_enterotype <- function(mapping = NULL, data = NULL, stat = "identity",  position = "identity", 
                            alpha = 0.3, prop = 0.5, ..., lineend = "butt", linejoin = "round", 
                            linemitre = 1, arrow = NULL, na.rm = FALSE, parse = FALSE, 
                            nudge_x = 0, nudge_y = 0, label.padding = unit(0.15, "lines"), 
                            label.r = unit(0.15, "lines"), label.size = 0.1, 
                            show.legend = TRUE, inherit.aes = TRUE) {
  library(ggplot2)
  # create new stat and geom for PCA scatterplot with ellipses
  StatEllipse <- ggproto("StatEllipse", Stat, 
                         required_aes = c("x", "y"), 
                         compute_group = function(., data, scales, level = 0.75, segments = 51, ...) {
                           library(MASS)
                           dfn <- 2
                           dfd <- length(data$x) - 1
                           if (dfd < 3) {
                             ellipse <- rbind(c(NA, NA))
                           } else {
                             v <- cov.trob(cbind(data$x, data$y))
                             shape <- v$cov
                             center <- v$center
                             radius <- sqrt(dfn * qf(level, dfn, dfd))
                             angles <- (0:segments) * 2 * pi/segments
                             unit.circle <- cbind(cos(angles), sin(angles))
                             ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
                           }
                           ellipse <- as.data.frame(ellipse)
                           colnames(ellipse) <- c("x", "y")
                           return(ellipse)
                         })
  
  # write new ggproto 
  GeomEllipse <- ggproto("GeomEllipse", Geom, 
                         draw_group = function(data, panel_scales, coord) {
                           n <- nrow(data)
                           if (n == 1) 
                             return(zeroGrob())
                           munched <- coord_munch(coord, data, panel_scales)
                           munched <- munched[order(munched$group), ]
                           first_idx <- !duplicated(munched$group)
                           first_rows <- munched[first_idx, ]
                           grid::pathGrob(munched$x, munched$y, default.units = "native", 
                                          id = munched$group, 
                                          gp = grid::gpar(col = first_rows$colour, 
                                                          fill = alpha(first_rows$fill, first_rows$alpha), lwd = first_rows$size * .pt, lty = first_rows$linetype))
                         }, 
                         default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1, alpha = NA, prop = 0.5), 
                         handle_na = function(data, params) {
                           data
                         }, 
                         required_aes = c("x", "y"), 
                         draw_key = draw_key_path
  )
  
  # create a new stat for PCA scatterplot with lines which totally directs to the center
  StatConline <- ggproto("StatConline", Stat, 
                         compute_group = function(data, scales) {
                           library(miscTools)
                           library(MASS)
                           df <- data.frame(data$x,data$y)
                           mat <- as.matrix(df)
                           center <- cov.trob(df)$center
                           names(center)<- NULL 
                           mat_insert <- insertRow(mat, 2, center )
                           for(i in 1:nrow(mat)) {
                             mat_insert <- insertRow( mat_insert, 2*i, center )
                             next
                           }
                           mat_insert <- mat_insert[-c(2:3),]
                           rownames(mat_insert) <- NULL
                           mat_insert <- as.data.frame(mat_insert,center)
                           colnames(mat_insert) =c("x","y")
                           return(mat_insert)
                         },
                         required_aes = c("x", "y")
                         
  )
  
  # create a new stat for PCA scatterplot with center labels
  StatLabel <- ggproto("StatLabel" ,Stat,
                       compute_group = function(data, scales) {
                         library(MASS)
                         df <- data.frame(data$x,data$y)
                         center <- cov.trob(df)$center
                         names(center)<- NULL 
                         center <- t(as.data.frame(center))
                         center <- as.data.frame(cbind(center))
                         colnames(center) <- c("x","y")
                         rownames(center) <- NULL
                         return(center)
                       },
                       required_aes = c("x", "y")
  )
  
  
  layer1 <- layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(na.rm = na.rm, ...))
  layer2 <- layer(stat = StatEllipse, data = data, mapping = mapping, geom = GeomEllipse, position = position, show.legend = FALSE, 
                  inherit.aes = inherit.aes, params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...))
  layer3 <- layer(data = data, mapping = mapping, stat =  StatConline, geom = GeomPath, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(lineend = lineend, linejoin = linejoin, 
                                linemitre = linemitre, arrow = arrow, na.rm = na.rm, ...))
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`", 
           call. = FALSE)
    }
    position <- position_nudge(nudge_x, nudge_y)
  }
  layer4 <- layer(data = data, mapping = mapping, stat = StatLabel, geom = GeomLabel, 
                  position = position, show.legend = FALSE, inherit.aes = inherit.aes, 
                  params = list(parse = parse, label.padding = label.padding, 
                                label.r = label.r, label.size = label.size, na.rm = na.rm, ...))
  return(list(layer1,layer2,layer3,layer4))
}

#source("individual_pca_functions.R")
# for subsets n=509
# vcf
setwd("~/Dropbox/Mac/Documents/HG/Domestication/01_pca")
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --remove UMFS_2_outlier.txt --recode --recode-INFO-all --out genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe", sep=""))

vcf.fn <- "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_509.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_509_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

n <- 17
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
#col_gradient = viridis_pal(option = "C")(17)  # n = number of colors seeked
#col_gradient = rep(c("#049DD9", "#F25C05"), c(9, 8))
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
tab_pop$pop <-factor(tab_pop$pop, levels=order1)
# for 2 clusters
order2 = c("DOM", "WILD")
tab_pop$pop1 <-factor(tab_pop$pop1, levels=order2)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop, shape=pop1), alpha = 0.8)+
  scale_color_manual(values=col_gradient , name="Population/Line") +
  scale_shape_manual(values=c(15, 17), name="Origin") +
  guides(fill = guide_legend(override.aes=list(shape=17)))+
  #scale_color_manual(values = c("#F25C05", "#049DD9"))+
  #scale_shape_manual(values=c(rep(c(0:2, 5, 15:18), 2), 6), breaks=order1)+
  #scale_shape_manual(values=c(rep(15, 9), rep(17,8)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+
  theme(legend.position="right",
        legend.title = element_text(size = 12),
        legend.text=element_text(size=12)) + 
  theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = "Black"))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

p2 <- p1 + stat_ellipse(aes(fill = pop), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  #scale_fill_manual(values = c("#F25C05", "#049DD9"))
  scale_fill_manual(values = col_gradient) +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(override.aes = list(shape = c(0,2)),
                              order = 1))
p2
graph2ppt(file="PCA_509",width=8,height=6)

jpeg("PCA_509.jpg", width = 8, height = 6, units = 'in', res = 300)
p2
dev.off()

#############################################################
# below is the PCA result without MEW and MEH1, total n=416 #
#############################################################

setwd("~/Dropbox/Mac/Documents/HG/Domestication/01_pca")
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --remove MEW_MEH1.txt --recode --recode-INFO-all --out genetyped_data_n_416_maf05_maxmiss095_popmiss095_hwe", sep=""))

vcf.fn <- "genetyped_data_n_416_maf05_maxmiss095_popmiss095_hwe.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_416_maf05_maxmiss095_popmiss095_hwe.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_416_maf05_maxmiss095_popmiss095_hwe.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_416_maf05_maxmiss095_popmiss095_hwe.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_416.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_416_edit.txt", header = TRUE, sep='\t')
tab_pop$population_nonum = stringr::str_remove(tab_pop$pop1, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop1) 

n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_gradient = viridis_pal(option = "C")(14)  # n = number of colors seeked
col_gradient = rep(c("#049DD9", "#F25C05"), c(7, 7))
col_gradient = rep(c("#049DD9", "#F25C05"), c(1, 1))
# PC1-2 for individual populations
order1 = c("LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2")
tab_pop$pop1 <-factor(tab_pop$pop, levels=order1)
# for 2 clusters
order1 = c("Dom", "Wild")
tab_pop$pop1 <-factor(tab_pop$pop1)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=4, aes(color=pop1, shape=pop1), alpha = 0.8)+
  labs(shape="Pop", colour="Pop")+
  scale_color_manual(values=col_gradient, breaks=order1)+
  #scale_color_manual(values = c("#F25C05", "#049DD9"))+
  #scale_shape_manual(values=c(rep(c(0:2, 5, 15:18), 2), 6), breaks=order1)+
  scale_shape_manual(values=c(rep(15, 9), rep(16,8)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text=element_text(size=12),
        text = element_text(size=14,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

p2 <- p1 + stat_ellipse(aes(fill = pop1), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  #scale_fill_manual(values = c("#F25C05", "#049DD9"))
  scale_fill_manual(values = col_gradient)
p2
graph2ppt(file="PCA_507_2",width=10,height=8)

jpeg("PCA_507.jpg", width = 10, height = 8, units = 'in', res = 300)
p2
dev.off()




