#############################################
# Check the af density of delta_p (wild)    #
#############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/Starting_allele_frequency_delta_p")

library(export)
# load reference file with header in it
ref = 'COH_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
ref_n = dat_ref$nInd
ref_k = round(dat_ref$knownEM*dat_ref$nInd*2)
# load challenge file with header in it
chr = 'ARN_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
ch_n <- dat_ch$nInd
ch_k <- round(dat_ch$knownEM*dat_ch$nInd*2)
# calculate delta_p
delta_p <- dat_ch$knownEM - dat_ref$knownEM
delta_p
length(delta_p)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
deltap_pos <- delta_p[delta_p > 0]
deltap_neg <- abs(delta_p[delta_p < 0])
# combine two datasets
dp <- c(deltap_neg, deltap_pos)
Directionality <- rep(c("negative", "positive"), times=c(length(deltap_neg),length(deltap_pos)))
df <- data.frame(dp,Directionality)
# start to plot
p1 <- df %>%
  ggplot( aes(x=dp, fill=Directionality)) +
  geom_histogram( color="#e9ecef", alpha=0.3, position = 'identity', bins = 40) +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ylim(0, 1.5e+05)+
  labs(fill = expression(Delta~italic("p")~directionality))+
  theme_classic()+
  ylab("Density") +
  xlab(expression(Delta~italic("p")~'('~ARN~'-'~COH~')'))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=20))+
  labs(fill="")
p1


ref = 'NB_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
ref_n = dat_ref$nInd
ref_k = round(dat_ref$knownEM*dat_ref$nInd*2)
# load challenge file with header in it
chr = 'HC_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
ch_n <- dat_ch$nInd
ch_k <- round(dat_ch$knownEM*dat_ch$nInd*2)
# calculate delta_p
delta_p <- dat_ch$knownEM - dat_ref$knownEM
delta_p
length(delta_p)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
deltap_pos <- delta_p[delta_p > 0]
deltap_neg <- abs(delta_p[delta_p < 0])
# combine two datasets
dp <- c(deltap_neg, deltap_pos)
Directionality <- rep(c("negative", "positive"), times=c(length(deltap_neg),length(deltap_pos)))
df <- data.frame(dp,Directionality)
# start to plot
p2 <- df %>%
  ggplot( aes(x=dp, fill=Directionality)) +
  geom_histogram( color="#e9ecef", alpha=0.3, position = 'identity', bins = 40) +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ylim(0, 1.5e+05)+
  labs(fill = expression(Delta~italic("p")~directionality))+
  theme_classic()+
  ylab("Density") +
  xlab(expression(Delta~italic("p")~'('~HC~'-'~NB~')'))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=20))
p2

ggarrange(p2, p1, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")

graph2ppt(file="ARN_COH_dp_density.pptx", width=8, height=6)

#########################################################
# Check the af distributin of starting allele (wild)    #
#########################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/Starting_allele_frequency_delta_p")
library(export)
# load reference file with header in it
ref = 'COH_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
# load challenge file with header in it
chr = 'ARN_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
# calculate delta_p
delta_p <- dat_ch$knownEM - dat_ref$knownEM
delta_p
length(delta_p)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
# calculate the number of SNPs at each interval when delta_p is negative
ref_af_neg <- dat_ref$knownEM[which(delta_p < 0)]
ref_af_neg <- ref_af_neg*100
neg <- as.data.frame(table(cut(ref_af_neg,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(neg)
# calculate the number of SNPs at each interval when delta_p is positive
ref_af_pos <- dat_ref$knownEM[which(delta_p > 0)]
ref_af_pos <- ref_af_pos*100
pos <- as.data.frame(table(cut(ref_af_pos,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(pos)
# combine two datasets
dat <- rbind(neg, pos)
dat$Directionality <- rep(c("negative", "positive"), each = 10)
# start to plot
p3 <- ggplot(dat, aes(fill=Directionality, y=Freq, x=Var1)) + 
  geom_bar(position="dodge", stat="identity")+
  ylim(0, 4e+05)+
  theme_classic()+
  ylab("Density") +
  xlab("Starting allele frequency (%) in COH")+
  labs(fill = expression(Delta~italic("p")~directionality))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1)) 
p3

ref = 'NB_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
# load challenge file with header in it
chr = 'HC_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
# calculate delta_p
delta_p <- dat_ch$knownEM - dat_ref$knownEM
delta_p
length(delta_p)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
# calculate the number of SNPs at each interval when delta_p is negative
ref_af_neg <- dat_ref$knownEM[which(delta_p < 0)]
ref_af_neg <- ref_af_neg*100
neg <- as.data.frame(table(cut(ref_af_neg,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(neg)
# calculate the number of SNPs at each interval when delta_p is positive
ref_af_pos <- dat_ref$knownEM[which(delta_p > 0)]
ref_af_pos <- ref_af_pos*100
pos <- as.data.frame(table(cut(ref_af_pos,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(pos)
# combine two datasets
dat <- rbind(neg, pos)
dat$Directionality <- rep(c("negative", "positive"), each = 10)
# start to plot
p4 <- ggplot(dat, aes(fill=Directionality, y=Freq, x=Var1)) + 
  geom_bar(position="dodge", stat="identity")+
  ylim(0, 4e+05)+
  theme_classic()+
  ylab("Density") +
  xlab("Starting allele frequency (%) in NB")+
  labs(fill = expression(Delta~italic("p")~directionality))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1)) 
p4
# convert to ppt format

ggarrange(p4, p3, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")

graph2ppt(file="ARN_COH_delta_p_dis_plot.pptx", width=8, height=6)

#########################################################
# Check the relationship between deltap and SGS p-value #
#########################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/Relationship_deltap_ps")
library(export)
# load reference file with header in it
dname = 'REF19_CHR19_ARN_COH_out_all_fish.txt'
dat <- read.delim(dname, header = TRUE, sep='\t')
delta_p <- dat$dp_w
length(dat$dp_w)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
# calculate the number of ps at each interval when delta_p is negative
dat_ps_neg <- dat$p_w[which(delta_p < 0)]
ps_neg <- dat_ps_neg*100
neg <- as.data.frame(table(cut(ps_neg,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(neg)
# calculate the number of ps at each interval when delta_p is positive
dat_ps_pos <- dat$p_w[which(delta_p > 0)]
ps_pos <- dat_ps_pos*100
pos <- as.data.frame(table(cut(ps_pos,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(pos)
# combine two datasets
dat <- rbind(neg, pos)
dat$Directionality <- rep(c("negative", "positive"), each = 10)
# start to plot
p5 <- ggplot(dat, aes(fill=Directionality, y=Freq, x=Var1)) + 
  geom_bar(position="dodge", stat="identity")+
  ylim(0, 250000)+
  theme_classic()+
  ylab("Density") +
  labs(fill = expression(Delta~italic("p")~directionality))+
  xlab(expression(SGS~italic("p")~values~'('~ARN~'-'~COH~')'))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1)) 
p5

dname = 'REF19_CHR19_NB_HC_out_all_fish.txt'
dat <- read.delim(dname, header = TRUE, sep='\t')
delta_p <- dat$dp_w
length(dat$dp_w)
length(delta_p[delta_p < 0])
# [1] 845337
length(delta_p[delta_p > 0])
# [1] 879182
# calculate the number of ps at each interval when delta_p is negative
dat_ps_neg <- dat$p_w[which(delta_p < 0)]
ps_neg <- dat_ps_neg*100
neg <- as.data.frame(table(cut(ps_neg,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(neg)
# calculate the number of ps at each interval when delta_p is positive
dat_ps_pos <- dat$p_w[which(delta_p > 0)]
ps_pos <- dat_ps_pos*100
pos <- as.data.frame(table(cut(ps_pos,seq(min(0),max(100),10),include.lowest = TRUE)))
plot(pos)
# combine two datasets
dat <- rbind(neg, pos)
dat$Directionality <- rep(c("negative", "positive"), each = 10)
# start to plot
p6 <- ggplot(dat, aes(fill=Directionality, y=Freq, x=Var1)) + 
  geom_bar(position="dodge", stat="identity")+
  ylim(0, 250000)+
  theme_classic()+
  ylab("Density") +
  labs(fill = expression(Delta~italic("p")~directionality))+
  xlab(expression(SGS~italic("p")~values~'('~HC~'-'~NB~')'))+
  theme(text = element_text(size=20))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1)) 
p6

ggarrange(p6, p5, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")

# convert to ppt format
graph2ppt(file="ARN_COH_ps_dis_plot.pptx", width=8, height=6)
