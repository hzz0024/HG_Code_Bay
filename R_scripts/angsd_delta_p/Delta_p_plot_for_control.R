####################################
##########  Delta_p plot ###########
####################################
setwd("~/Documents/Ryan_workplace/DelBay_adult/08_fish_exact")
source("manhattan.R")
###################### load maf values for delta_p ##########

file1 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_REF20_CHR19_CHR20_out_0.05_696.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_control.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Ref - 2020 Ref '*Delta~italic(p)), cex.lab=1.5) 

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'CHR20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_REF20_CHR19_CHR20_out_0.05_696.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

#jpeg("Mahattan_ch5_dxy_150.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(dat1, chr="chromo",bp="position",p="deltaP", snp = "SNP", highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2020 Sur '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, col = "grey80")
dev.off()

###################### test on chr1 ##########
source("manhattan.R")
file1 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_REF20_CHR19_CHR20_outlier_chr1_16.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_control_chr1.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 1), xlim = c(8000000, 10000000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Ref - 2020 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'CHR20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
#file2 = 'REF19_REF20_CHR19_CHR20_out_0.05_696.csv' 
file2 = 'REF19_REF20_CHR19_CHR20_outlier_chr1_16.csv'
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 1), xlim = c(8000000, 10000000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2020 Sur '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

dev.off()


###################### test on chr5 ##########
source("manhattan.R")
file1 = 'REF19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'REF20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
file2 = 'REF19_REF20_CHR19_CHR20_outlier_chr5_26.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

jpeg("Mahattan_control_chr5.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(10000000, 16000000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Ref - 2020 Ref '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")

file1 = 'CHR19_minq20_minmq30_1x_CV30_masked.mafs' #p1 CH
file0 = 'CHR20_minq20_minmq30_1x_CV30_masked.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM 
p0 = dat0$knownEM 
delta_p = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p
#file2 = 'REF19_REF20_CHR19_CHR20_out_0.05_696.csv' 
file2 = 'REF19_REF20_CHR19_CHR20_outlier_chr5_26.csv'
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(10000000, 16000000), highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression('2019 Sur - 2020 Sur '*Delta~italic(p)), cex.lab=1.5) 
abline(h=0, lty = 2, col = "red")
dev.off()