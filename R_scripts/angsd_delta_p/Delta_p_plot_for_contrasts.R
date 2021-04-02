####################################
##########  Delta_p plot ###########
####################################
setwd("~/Documents/Ryan_workplace/DelBay_adult/08_fish_exact")
source("manhattan.R")
###################### load maf values for delta_p ##########

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_CHR19_NB_HC_Z_outlier_FDR0.05.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_REF19_CHR19_NB_HC.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2019 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'NB_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_CHR19_NB_HC_Z_outlier_FDR0.05.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - NB '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()


###################### REF20_CHR20_NB_HC_Z_outlier_FDR0.05 ###############

###################### load maf values for delta_p ##########

file1 = 'CHR20_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF20_CHR20_NB_HC_Z_outlier_FDR0.05.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_REF20_CHR20_NB_HC.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2020 Sur - 2019 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'NB_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF20_CHR20_NB_HC_Z_outlier_FDR0.05.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - NB '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()


###################### NC_035784.1_16552716 ##########

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_NC_035784.1_16552716.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2019 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'NB_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - NB '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'SR_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - SR '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()



###################### SGS outlier delta_p ##########

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'SGS_outlier_3006.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_SGS.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(1,1))
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2019 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()



###################### Mahattan_PRF.jpg ##########
source("manhattan.R")

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
file1 = 'NC_035784.1_16552962.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)

jpeg("Mahattan_PRF.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2019 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'NB_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
file1 = 'NC_035784.1_16552962.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)

manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight1 = h1$V1,highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - NB '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'HC_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'SR_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'NC_035784.1_16552716.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)
file1 = 'NC_035784.1_16552962.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)

#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16500000, 16600000), highlight1 = h1$V1,highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 HC - SR '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()