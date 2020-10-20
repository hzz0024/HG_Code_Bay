library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)
# load the R function
source("manhattan.R")
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
###################### 150 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest2.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(16534778, 16560000),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1, main = "CH-REF chr5:150 SNP/window ",) 

###################### start to add Tajima's D ####################

# REF theta and tajima'D
t_name1 = 'REF.1000_thetas.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
# CH theta and tajima'D
t_name2 = 'CH.1000_thetas.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')

idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
yv1 = t1$Tajima[idx1]
low1 = floor(min(yv1))
high1 = ceiling((max(yv1)))

idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
yv2 = t2$Tajima[idx2]
low2 = floor(min(yv2))
high2 = ceiling((max(yv2)))

low = min(c(low1,low2))
high = max(c(high1, high2))

par(mar=c(5,4,4,5)+.1)
axis(4, at=seq(0.0,0.5,0.1),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=0.8)
mtext('Tajima\'s D', side=4, line=3, cex=1)

lines(xv1, (yv1-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')
points(xv1,(yv1-low)/(high-low)*(0.5-0), col='blue')
lines(xv2, (yv2-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='orange')
points(xv2,(yv2-low)/(high-low)*(0.5-0), col='orange')

legend(16553000, 0.5, legend=c("REF Tajima's D", "CH Tajima's D"), col=c("blue", "orange"), lty=1, cex=0.8)


###################### start new delta_p plot ###############
###################### load maf values for delta_p ##########
source("manhattan.R")
file1 = 'CH_maf0.05_pctind0.7_cv30.mafs' #p1 CH
file0 = 'REF_maf0.05_pctind0.7_cv30.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
###################### 150 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest2.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
file2 = 'outliers.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16534778, 16560000), highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1, main = "CH-REF chr5:150 SNP/window ",) 
# draw lines
#idx = dat1$chromo == 5
#lines(dat1$position[idx], delta_p[idx],col='gray90')
###################### start to add Tajima's D ####################

# REF theta and tajima'D
t_name1 = 'REF.1000_thetas.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
# CH theta and tajima'D
t_name2 = 'CH.1000_thetas.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')

idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
yv1 = t1$Tajima[idx1]
low1 = floor(min(yv1))
high1 = ceiling((max(yv1)))

idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
yv2 = t2$Tajima[idx2]
low2 = floor(min(yv2))
high2 = ceiling((max(yv2)))

low = min(c(low1,low2))
high = max(c(high1, high2))

par(mar=c(5,4,4,5)+.1)
axis(4, at=seq(-0.5,0.5,0.2),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=0.8)
mtext('Tajima\'s D', side=4, line=3, cex=1)

lines(xv1, (yv1-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')
points(xv1,(yv1-low)/(high-low)*(0.5-0), col='blue')
lines(xv2, (yv2-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='orange')
points(xv2,(yv2-low)/(high-low)*(0.5-0), col='orange')

legend(16553000, 0.5, legend=c("REF Tajima's D", "CH Tajima's D"), col=c("blue", "orange"), lty=1, cex=0.8)

########################### start theta plot ##########################



###################### 150 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest2.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(DT1, chromo == 5), xlim = c(16534778, 16560000),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF dxy", cex.lab=1, main = "CH-REF chr5:150 SNP/window ",) 

###################### start to add theta ####################

# REF theta and tajima'D
t_name1 = 'REF.1000_thetas.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
# CH theta and tajima'D
t_name2 = 'CH.1000_thetas.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')

idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
theta1 = t1$tP/t1$nSites
theta1[is.na(theta1)] <- 0
#yv1 = t1$tW[idx1]
yv1 = theta1[idx1]
low1 = floor(min(yv1))
high1 = ceiling((max(yv1)))

idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
theta2 = t2$tP/t2$nSites
theta2[is.na(theta2)] <- 0
#yv2 = t2$tW[idx2]
yv2 = theta2[idx2]
low2 = floor(min(yv2))
high2 = ceiling((max(yv2)))

low = min(c(low1,low2))
high = max(c(high1, high2))

#low = 0
#high = 0.2

par(mar=c(5,4,4,5)+.1)
axis(4, at=seq(0.0,0.5,0.1),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=0.8)
mtext('nucleotide diversity', side=4, line=3, cex=1)

lines(xv1, (yv1-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')
#points(xv1,(yv1-low)/(high-low)*(0.5-0), col='blue')
lines(xv2, (yv2-low)/(high-low)*(0.5-0), xaxt='n',yaxt='n', xlab='', ylab='', col='orange')
#points(xv2,(yv2-low)/(high-low)*(0.5-0), col='orange')

legend(16554000, 0.5, legend=c("REF theta", "CH theta"), col=c("blue", "orange"), lty=1, cex=0.8)


