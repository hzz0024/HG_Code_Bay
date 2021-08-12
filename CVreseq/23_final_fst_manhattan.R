source("manhattan.R")
library(KRIS)
options(scipen=999)
library(export)

headname = "CS_HC_noinvers." 
titlename = "CS_HC"

win = 1000
name = paste0(headname, win, "bp.", "s", win/5, ".csv")
DT = read.delim(name, header = TRUE, sep=',')
mid_pos <- round((DT$start + DT$end)/2)
id = paste0(DT$scaffold,'_',mid_pos)
DT <- as.data.frame(cbind(DT,mid_pos, id))
DT <- DT[complete.cases(DT), ]
DT[,9][DT[,9]<0] = 0.000001 #@chnage
zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE))# @change
dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9], zfst = zfst) # @change
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$fst <- as.numeric(dat$fst)
dat$zfst <- as.numeric(dat$zfst)
# calculate the 99.9% quantile
thred=quantile(dat$zfst, 0.999)
cnt = length(dat[dat$zfst>thred[[1]],1])
outlier = dat[dat$zfst>thred[[1]],5]
#cnt = sum(dat$zfst>5)
#outlier = dat[dat$zfst>5,5]
print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))


file1 = 'CS_HC_fst_snpofinterest.txt' 
h1 = read.delim(file1, header = FALSE, sep='\t')
h1 = as.list(h1)

jpeg("manhattan1.jpg", width = 12, height = 6, units = 'in', res = 300)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(2,1))
manhattan(dat, chr="chr",bp="mid_pos",p="zfst", highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, max(dat$zfst)+0.2), #subset(dat, chr == 8)
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression("Z"~italic(F)[ST]), cex.lab=1.5, xlab = NA, cex.main=1.5)
#graph2ppt(file="manhattan1",width=6,height=10)
abline(h=thred, col='red', lwd=1)
dev.off()

headname = "CLP_HCVA_noinvers." 
titlename = "CLP_HCVA"

win = 1000
name = paste0(headname, win, "bp.", "s", win/5, ".csv")
DT = read.delim(name, header = TRUE, sep=',')
mid_pos <- round((DT$start + DT$end)/2)
id = paste0(DT$scaffold,'_',mid_pos)
DT <- as.data.frame(cbind(DT,mid_pos, id))
DT <- DT[complete.cases(DT), ]
DT[,9][DT[,9]<0] = 0.000001 #@chnage
zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE))# @change
dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9], zfst = zfst) # @change
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$fst <- as.numeric(dat$fst)
dat$zfst <- as.numeric(dat$zfst)
# calculate the 99.9% quantile
thred=quantile(dat$zfst, 0.999)
cnt = length(dat[dat$zfst>thred[[1]],1])
outlier = dat[dat$zfst>thred[[1]],5]
#cnt = sum(dat$zfst>5)
#outlier = dat[dat$zfst>5,5]
print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))

file1 = 'CLP_HCVA_fst_snpofinterest.txt' 
h1 = read.delim(file1, header = FALSE, sep='\t')
h1 = as.list(h1)

jpeg("manhattan2.jpg", width = 12, height = 6, units = 'in', res = 300)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(2,1))
manhattan(dat, chr="chr",bp="mid_pos",p="zfst", highlight2 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, max(dat$zfst)+0.2), #subset(dat, chr == 8)
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab=expression("Z"~italic(F)[ST]), cex.lab=1.5, cex.main=1.5)
abline(h=thred, col='red', lwd=1)
#graph2ppt(file="manhattan1",width=6,height=10)
dev.off()

