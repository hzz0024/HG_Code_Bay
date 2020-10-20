library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)
# load the R function
source("manhattan.R")
# load snp of interest
file1 = 'snpsOfInterest.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
file2 = 'outliers.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
# load the CH-REF dxy
name1 = "CH_REF_1934038_dxy"
DT1 = read.delim(name1, header = TRUE, sep='\t')
id1 = paste0(DT1$chromo,'_',DT1$position)
DT1$SNP <- id1 
# load the HC_SR dxy
name2 = "HC_SR_1934038_dxy"
DT2 = read.delim(name2, header = TRUE, sep='\t')
id2 = paste0(DT2$chromo,'_',DT2$position)
DT2$SNP <- id2
# load the HC_NB dxy
name3 = "HC_NB_1934038_dxy"
DT3 = read.delim(name3, header = TRUE, sep='\t')
id3 = paste0(DT3$chromo,'_',DT3$position)
DT3$SNP <- id3
##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)
DT1$chromo <- as.numeric(DT1$chromo)
DT2$chromo <- as.numeric(DT2$chromo)
DT3$chromo <- as.numeric(DT3$chromo)
###################### 10K SNPs/window ######################
jpeg("Mahattan_ch5_dxy_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
                    col=c("grey","black"),genomewideline=F, suggestiveline=F,
                    ylab="CH_REF dxy", cex.lab=1.4, main = "CH-REF chr5:10K SNP/window ",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT2, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
                    col=c("grey","black"),genomewideline=F, suggestiveline=F,
                    ylab="HC-SR dxy", cex.lab=1.4, main = "HC-SR chr5:10K SNP/window",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT3, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
                    col=c("grey","black"),genomewideline=F, suggestiveline=F,
                    ylab="HC-NB dxy", cex.lab=1.4, main = "HC-NB chr5:10K SNP/window",) 
dev.off()
###################### 1500 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest4.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1.4, main = "CH-REF chr5:1500 SNP/window ",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT2, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR dxy", cex.lab=1.4, main = "HC-SR chr5:1500 SNP/window",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT3, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB dxy", cex.lab=1.4, main = "HC-NB chr5:1500 SNP/window",) 
dev.off()
###################### 150 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest2.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1.4, main = "CH-REF chr5:150 SNP/window ",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT2, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR dxy", cex.lab=1.4, main = "HC-SR chr5:150 SNP/window",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT3, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB dxy", cex.lab=1.4, main = "HC-NB chr5:150 SNP/window",) 
dev.off()
###################### 25 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest3.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_dxy_25.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1.4, main = "CH-REF chr5:25 SNP/window ",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT2, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR dxy", cex.lab=1.4, main = "HC-SR chr5:25 SNP/window",) 
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT3, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB dxy", cex.lab=1.4, main = "HC-NB chr5:25 SNP/window",) 
dev.off()



