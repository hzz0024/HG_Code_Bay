library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)
# load the R function

# set up the working directory
setwd("~/Documents/HG/DelBay19_adult/09_theta/plot")
# source the manhattan plot function
source("manhattan.R")

########## process delta_p ########## 

file1 = 'CHR19_all_minq20_minmq30_CV30_masked.mafs' #p1 CH
file0 = 'REF19_all_minq20_minmq30_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
dat <- cbind(dat1$chromo, dat1$position, dat1$deltaP)
write.table(dat, file = "CHR19_REF19_delta_p.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

########## manhattan ##########
source("manhattan.R")
# load snp of interest
file1 = '5_16552716.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_TajimaD_diff_challenge_ch5_5_16552716.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(7,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('Surv - Ref '*Delta~italic(p)), cex.lab=1.4) 

#manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 1), xlim = c(1, 100000), highlight2 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
#          col=c("grey","black"),genomewideline=F, suggestiveline=F,
#          ylab=expression('Surv - Ref '*Delta~italic(p)), cex.lab=1.4) 

abline(h=0, col="grey", lty=2)
###################### start to add Tajima's D ####################

# REF theta and tajima'D
t_name1 = 'REF19.thetas.1000.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
#t1 = t1[complete.cases(t1), ]
# extract SNP from chr5
idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
# extract ranges
idx=xv1>16500000 & xv1<16600000
xv1 = xv1[idx]
yv1 = t1$Tajima[idx1]
yv1 = yv1[idx]

# CH theta and tajima'D
t_name2 = 'CHR19.thetas.1000.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')
#t2 = t2[complete.cases(t2), ]
low1 = floor(min(yv1))
high1 = ceiling((max(yv1)))
idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
# extract ranges
idx=xv2>16500000 & xv2<16600000
xv2 = xv2[idx]
yv2 = t2$Tajima[idx2]
yv2 = yv2[idx]

tajima_diff = yv2 - yv1
dat = cbind(xv1, xv2, tajima_diff)
dat = dat[complete.cases(dat), ]

low = floor(min(dat[,3]))
high = ceiling((max(dat[,3])))

#high = as.numeric(format(round(high, 2), nsmall = 2))

par(mar=c(7,8,1,5)+.1)
#par(mar=c(4,4,4,5)+.1)
axis(4, at=seq(-0.5,0.5,0.2),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=1.2)
mtext('Tajima\'s D (Surv - Ref)', side=4, line=3, cex=1.2)
# start scaling
lines(dat[,1], (dat[,3]-low)/(high-low)*(0.5-(-0.5))+(-0.5), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')

#legend(16585000, -0.45, legend=c("Reference Tajima's D", "Survivors Tajima's D"), col=c("blue", "orange"), lty=1, cex=0.8)
dev.off()


###################### second plot for deltap and diversity ####################
jpeg("Mahattan_pi_diff_challenge_ch5_5_16552716.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(7,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('Surv - Ref '*Delta~italic(p)), cex.lab=1.4) 
abline(h=0, col="grey", lty=2)

###################### start to add theta ####################
source("manhattan.R")
# load snp of interest
file1 = '5_16552716.csv' 
h2 = read.delim(file1, header = FALSE, sep=',')
h2 = as.list(h1)
# REF theta and tajima'D
t_name1 = 'REF19.thetas.1000.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
#t1 = t1[complete.cases(t1), ]

idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
idx=xv1>16500000 & xv1<16600000
xv1 = xv1[idx]
theta1 = t1$tP/t1$nSites
theta1[is.na(theta1)] <- 0
#yv1 = t1$tW[idx1]
yv1 = theta1[idx1]
yv1 = yv1[idx]

# CH theta and tajima'D
t_name2 = 'CHR19.thetas.1000.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')
#t2 = t2[complete.cases(t2), ]
idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
idx=xv2>16500000 & xv2<16600000
xv2 = xv2[idx]
theta2 = t2$tP/t2$nSites
theta2[is.na(theta2)] <- 0
#yv2 = t2$tW[idx2]
yv2 = theta2[idx2]
yv2 = yv2[idx]

pi_diff = yv2 - yv1
dat = cbind(xv1, xv2, pi_diff)
dat = dat[complete.cases(dat), ]

low = as.numeric(formatC(min(dat[,3]), digits = 3, format = "f"))
high = as.numeric(formatC(max(dat[,3]), digits = 3, format = "f"))

par(mar=c(7,8,1,5)+.1)
axis(4, at=seq(-0.5,0.5,0.2),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=0.8)
mtext('Nucleotide diversity (Surv - Ref)', side=4, line=3, cex=1)

lines(dat[,1], (dat[,3]-low)/(high-low)*(0.5-(-0.5))+(-0.5), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')
dev.off()

########## process wild populations ########## 

file1 = 'HC_all_minq20_minmq30_CV30_masked.mafs' #p1 CH
file0 = 'NB_all_minq20_minmq30_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #HC
p0 = dat0$knownEM #NB
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
dat <- cbind(dat1$chromo, dat1$position, dat1$deltaP)
write.table(dat, file = "HC_NB_delta_p.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

########## manhattan ##########
source("manhattan.R")
# load snp of interest
file1 = '5_16552716.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_TajimaD_diff_wild_ch5_5_16552716.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(7,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('HC - NB '*Delta~italic(p)), cex.lab=1.4) 

#manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 1), xlim = c(1, 100000), highlight2 = h1$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
#          col=c("grey","black"),genomewideline=F, suggestiveline=F,
#          ylab=expression('Surv - Ref '*Delta~italic(p)), cex.lab=1.4) 

abline(h=0, col="grey", lty=2)
###################### start to add Tajima's D ####################

# NB theta and tajima'D
t_name1 = 'NB.thetas.1000.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
#t1 = t1[complete.cases(t1), ]
# extract SNP from chr5
idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
# extract ranges
idx=xv1>16500000 & xv1<16600000
xv1 = xv1[idx]
yv1 = t1$Tajima[idx1]
yv1 = yv1[idx]

# HC theta and tajima'D
t_name2 = 'HC.thetas.1000.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')
#t2 = t2[complete.cases(t2), ]
low1 = floor(min(yv1))
high1 = ceiling((max(yv1)))
idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
# extract ranges
idx=xv2>16500000 & xv2<16600000
xv2 = xv2[idx]
yv2 = t2$Tajima[idx2]
yv2 = yv2[idx]

tajima_diff = yv2 - yv1
dat = cbind(xv1, xv2, tajima_diff)
dat = dat[complete.cases(dat), ]

low = floor(min(dat[,3]))
high = ceiling((max(dat[,3])))

#high = as.numeric(format(round(high, 2), nsmall = 2))

par(mar=c(7,8,1,5)+.1)
#par(mar=c(4,4,4,5)+.1)
axis(4, at=seq(-0.5,0.5,0.2),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=1.2)
mtext('Tajima\'s D (HC - NB)', side=4, line=3, cex=1.2)
# start scaling
lines(dat[,1], (dat[,3]-low)/(high-low)*(0.5-(-0.5))+(-0.5), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')

#legend(16585000, -0.45, legend=c("Reference Tajima's D", "Survivors Tajima's D"), col=c("blue", "orange"), lty=1, cex=0.8)
dev.off()


###################### second plot for deltap and diversity ####################
jpeg("Mahattan_theta_diff_wild_ch5_5_16552716.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(7,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('HC - NB '*Delta~italic(p)), cex.lab=1.4) 
abline(h=0, col="grey", lty=2)

###################### start to add theta ####################
source("manhattan.R")
# load snp of interest
file1 = '5_16552716.csv' 
h2 = read.delim(file1, header = FALSE, sep=',')
h2 = as.list(h1)
# REF theta and tajima'D
t_name1 = 'NB.thetas.1000.window.idx.pestPG' 
t1 = read.delim(t_name1, header = TRUE, sep='\t')
#t1 = t1[complete.cases(t1), ]

idx1 = t1$Chr == 5
xv1 = t1$WinCenter[idx1]
idx=xv1>16500000 & xv1<16600000
xv1 = xv1[idx]
theta1 = t1$tP/t1$nSites
theta1[is.na(theta1)] <- 0
#yv1 = t1$tW[idx1]
yv1 = theta1[idx1]
yv1 = yv1[idx]

# CH theta and tajima'D
t_name2 = 'HC.thetas.1000.window.idx.pestPG' 
t2 = read.delim(t_name2, header = TRUE, sep='\t')
#t2 = t2[complete.cases(t2), ]
idx2 = t2$Chr == 5
xv2 = t2$WinCenter[idx2]
idx=xv2>16500000 & xv2<16600000
xv2 = xv2[idx]
theta2 = t2$tP/t2$nSites
theta2[is.na(theta2)] <- 0
#yv2 = t2$tW[idx2]
yv2 = theta2[idx2]
yv2 = yv2[idx]

pi_diff = yv2 - yv1
dat = cbind(xv1, xv2, pi_diff)
dat = dat[complete.cases(dat), ]

low = as.numeric(formatC(min(dat[,3]), digits = 3, format = "f"))
high = as.numeric(formatC(max(dat[,3]), digits = 3, format = "f"))

par(mar=c(7,8,1,5)+.1)
axis(4, at=seq(-0.5,0.5,0.2),seq(low,high,(high-low)/5),col.ticks='black',col.axis='black', las=1, cex.axis=0.8)
mtext('Nucleotide diversity (HC - NB)', side=4, line=3, cex=1)

lines(dat[,1], (dat[,3]-low)/(high-low)*(0.5-(-0.5))+(-0.5), xaxt='n',yaxt='n', xlab='', ylab='', col='blue')
dev.off()



