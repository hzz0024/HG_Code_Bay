#############################################
#########  output the outlier id   ##########
#############################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/ngsld_theta_distribution")
pname = "ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
head(dat)
df1 <- dat[which(dat$adj< 0.05),1:2]
df2 <- paste0(df1$chromo, "\t", df1$position, "\t")
message(paste0("number of outlier in ", pname, " is ", length(df$chromo)))
write.table(df2, "./ps_Del19_HC_NB_outlier_FDR05.txt", row.names=F, col.names=F, quote=F, sep="\t")

################################################################
#########  plot the theta distribution for outliers   ##########
################################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/ngsld_theta_distribution")
pname = "test.output.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
colnames(dat) = c('win', 'Chr', 'WinCenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'Tajima', 'fuf', 'fud', 'fayh', 'zeng', 'nSites')
dat$theta <- dat$tP/dat$nSites
hist(dat$theta, probability = T)
lines(density(dat$theta), col=2)
# test for Normality
ks.test(dat$theta, "pnorm", mean=mean(dat$theta), sd=sd(dat$theta))
# find the 95% CI
summary(dat$theta)
t.test(dat$theta, conf.level = 0.95)

mean(dat$theta)
sd(dat$theta)

#########################################################################
#########  random sample snps from the genome with similar pi  ##########
#########################################################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/ngsld_theta_distribution")
pname = "sampled.NB.theta.window.idx.pestPG"
dat = read.delim(pname, header = FALSE, sep='\t')
colnames(dat) = c('win', 'Chr', 'WinCenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'Tajima', 'fuf', 'fud', 'fayh', 'zeng', 'nSites')
dat$theta <- dat$tP/dat$nSites
dat <- dat[complete.cases(dat), ]
hist(dat$theta)
# test for Normality
ks.test(dat$theta, "pnorm", mean=mean(dat$theta), sd=sd(dat$theta))
# find the 95% CI
summary(dat$theta)
t.test(dat$theta, conf.level = 0.95)

mean(dat$theta)
sd(dat$theta)
