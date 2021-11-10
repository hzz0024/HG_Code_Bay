#############################################
# A simple comp between angular vs original #
#############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS/11_manhattan/")
pname = "ps_Del19_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
outlier_cnt <- length(dat$adj[which(dat$adj< 0.05)])
outlier <- dat[which(dat$adj< 0.05),]
dat$FDR = -log(dat$adj, 2)
message("number of outliers is ", outlier_cnt)
id1 = paste0(outlier$chromo,'_',outlier$position)

source("manhattan.R")
jpeg("Mahattan_CHR19_REF19.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="FDR", snp = "SNP", dat, highlight1 = id1 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 6),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log2(FDR)", cex.lab=1.5) 
dev.off()

pname = "ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
outlier_cnt <- length(dat$adj[which(dat$adj< 0.05)])
outlier <- dat[which(dat$adj< 0.05),]
dat$FDR = -log(dat$adj, 2)
message("number of outliers is ", outlier_cnt)
id1 = paste0(outlier$chromo,'_',outlier$position)

source("manhattan.R")
jpeg("Mahattan_HC_NB.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="FDR", snp = "SNP", dat, highlight1 = id1 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 6),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log2(FDR)", cex.lab=1.5) 
dev.off()

pname = ""
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
outlier_cnt <- length(dat$adj[which(dat$adj< 0.05)])
outlier <- dat[which(dat$adj< 0.05),]
dat$FDR = -log(dat$adj, 2)
message("number of outliers is ", outlier_cnt)
id1 = paste0(outlier$chromo,'_',outlier$position)

source("manhattan.R")
jpeg("Mahattan_HC_NB.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(1,1))
manhattan(chr="chromo",bp="position",p="FDR", snp = "SNP", dat, highlight1 = id1 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 6),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log2(FDR)", cex.lab=1.5) 
dev.off()
