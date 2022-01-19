####################
# Examine the data #
####################
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

# see https://github.com/clairemerot/Intro2PopGenomics/blob/master/3.2.3/PCA/script_PCA_from_vcf.R for original code
library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

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
