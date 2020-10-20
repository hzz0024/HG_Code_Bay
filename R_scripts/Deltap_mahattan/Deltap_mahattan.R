source("manhattan.R")

# CH-REF
file1 = 'CH_maf0.05_pctind0.7_cv30.mafs' #p1 CH
file0 = 'REF_maf0.05_pctind0.7_cv30.mafs' #p0 REF
dat1 = read.delim(file1, header = TRUE, sep='\t')
dat0 = read.delim(file0, header = TRUE, sep='\t')
p1 = dat1$knownEM #CH
p0 = dat0$knownEM #REF
delta_p1 = p1-p0
id = paste0(dat0$chromo,'_',dat0$position)
dat1$SNP <- id
dat1$deltaP <- delta_p1

# HC-SR
file3 = 'HC_maf0.05_pctind0.7_cv30.mafs' #p3 HC
file2 = 'SR_maf0.05_pctind0.7_cv30.mafs' #p2 SR
dat3 = read.delim(file3, header = TRUE, sep='\t')
dat2 = read.delim(file2, header = TRUE, sep='\t')
p3 = dat3$knownEM #HC
p2 = dat2$knownEM #SR
delta_p2 = p3-p2
id = paste0(dat2$chromo,'_',dat2$position)
dat2$SNP <- id
dat2$deltaP <- delta_p2

# HC-NB
file5 = 'HC_maf0.05_pctind0.7_cv30.mafs' #p5 HC
file4 = 'NB_maf0.05_pctind0.7_cv30.mafs' #p4 NB
dat5 = read.delim(file5, header = TRUE, sep='\t')
dat4 = read.delim(file4, header = TRUE, sep='\t')
p5 = dat5$knownEM #HC
p4 = dat4$knownEM #SR
delta_p3 = p5-p4
id = paste0(dat3$chromo,'_',dat3$position)
dat3$SNP <- id
dat3$deltaP <- delta_p3

file2 = 'outliers.csv' 
h2 = read.delim(file2, header = FALSE, sep=',')
h2 = as.list(h2)

###################### 10K SNPs/window ######################
file1 = 'snpsOfInterest.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_delta_p_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF delta_p", cex.lab=1.4, main = "CH-REF chr5:10K SNP/window ",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat2, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR delta_p", cex.lab=1.4, main = "HC-SR chr5:10K SNP/window",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat3, chromo == 5), xlim = c(13888433, 18446624),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB delta_p", cex.lab=1.4, main = "HC-NB chr5:10K SNP/window",) 
abline(h=0, col = "grey80")
dev.off()
###################### 1500 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest4.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_delta_p_1500.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF delta_p", cex.lab=1.4, main = "CH-REF chr5:1500 SNP/window ",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat2, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR delta_p", cex.lab=1.4, main = "HC-SR chr5:1500 SNP/window",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat3, chromo == 5), xlim = c(16377719, 16669307),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB delta_p", cex.lab=1.4, main = "HC-NB chr5:1500 SNP/window",) 
abline(h=0, col = "grey80")
dev.off()
###################### 150 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest2.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_delta_p_150.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF delta_p", cex.lab=1.4, main = "CH-REF chr5:150 SNP/window ",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat2, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR delta_p", cex.lab=1.4, main = "HC-SR chr5:150 SNP/window",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat3, chromo == 5), xlim = c(16534778, 16563728),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB delta_p", cex.lab=1.4, main = "HC-NB chr5:150 SNP/window",) 
abline(h=0, col = "grey80")
dev.off()
###################### 25 SNPs/window ######################
# load snp of interest
file1 = 'snpsOfInterest3.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
jpeg("Mahattan_ch5_delta_p_25.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat1, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="CH_REF delta_p", cex.lab=1.4, main = "CH-REF chr5:25 SNP/window ",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat2, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-SR delta_p", cex.lab=1.4, main = "HC-SR chr5:25 SNP/window",) 
abline(h=0, col = "grey80")
manhattan(chr="chromo",bp="position",p="deltaP", snp = "SNP", subset(dat3, chromo == 5), xlim = c(16550548, 16554169),  highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 0.8, ylim = c(-0.6, 0.6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="HC-NB delta_p", cex.lab=1.4, main = "HC-NB chr5:25 SNP/window",) 
abline(h=0, col = "grey80")
dev.off()



