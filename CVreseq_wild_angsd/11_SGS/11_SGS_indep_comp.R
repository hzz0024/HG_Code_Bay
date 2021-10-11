#############################################
# A simple comp between angular vs original #
#############################################
setwd("~/Dropbox/Mac/Documents/HG/CVreseq_wild_angsd/11_SGS/11_comp_lcwg_wgs")
pname = "ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
#dat$V6[dat$V6 == 0] = 0.00001
#dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP')

cnt_outlier <- length(dat$ps[which(dat$ps< 0.01)])
angular_outlier <- dat[which(dat$ps< 0.01),]
id1 <- paste0(angular_outlier$chromo,'_',angular_outlier$position)

neg <- length(angular_outlier$ps[which(angular_outlier$delta_p<0)])
pos <-  length(angular_outlier$ps[which(angular_outlier$delta_p>=0)])
message("number of outliers with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

neg <- length(dat$delta_p[which(dat$delta_p<0)])
pos <- length(dat$delta_p[which(dat$delta_p>=0)])
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

dat$ps = -log(dat$ps)

source("manhattan.R")
jpeg("Mahattan_adj_comp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(chr="chromo",bp="position",p="ps", snp = "SNP", subset(dat, chromo == 10), highlight1 = id1 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 15),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log(p-value)", cex.lab=1.5) 

pname = "ps_HC_CS19.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
#dat$V6[dat$V6 == 0] = 0.00001
#dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP')

cnt_outlier <- length(dat$ps[which(dat$ps< 0.01)])
angular_outlier <- dat[which(dat$ps< 0.01),]
id2 <- paste0(angular_outlier$chromo,'_',angular_outlier$position)

neg <- length(angular_outlier$ps[which(angular_outlier$delta_p<0)])
pos <-  length(angular_outlier$ps[which(angular_outlier$delta_p>=0)])
message("number of outliers with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

neg <- length(dat$delta_p[which(dat$delta_p<0)])
pos <- length(dat$delta_p[which(dat$delta_p>=0)])
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

dat$ps = -log(dat$ps)

manhattan(chr="chromo",bp="position",p="ps", snp = "SNP", subset(dat, chromo == 10), highlight1 = id2 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 15),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log(p-value)", cex.lab=1.5) 

dev.off()

length(intersect(id1,id2))
# 0
# > length(id2) # WGS
# [1] 2268
# > length(id1) # lcwg
# [1] 2268
