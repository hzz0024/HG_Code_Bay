# https://stats.stackexchange.com/questions/92124/how-to-interpret-a-qq-plot-of-p-values
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS/11_QQplot_pvalues")

simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)), xlab="Expect p-value (-log10 scale)", ylab="Observed p-value (-log10 scale)")
  abline(0, 1, col = "red")
  axis(2,cex.axis=2)
}

pname = "ps_Del19_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
head(dat)
df1 <- dat[which(dat$adj< 0.05),]
message(paste0("number of outlier in ", pname, " is ", length(df1$chromo)))

par(mar=c(4,8,1,6))
jpeg("QQplot_adj.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
simpleQQPlot(dat$adj)
dev.off()

par(mar=c(4,8,1,6))
jpeg("QQplot_ps.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
simpleQQPlot(dat$ps)
dev.off()
