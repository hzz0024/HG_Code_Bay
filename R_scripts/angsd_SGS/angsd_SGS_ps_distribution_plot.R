library("ggplot2")
library("plyr")
library(gridExtra)
setwd("~/Documents/Ryan_workplace/DelBay_adult/11_SGS")
# load the ps datasets
dat1 = read.delim('ps_Del19_challenge.txt', header = FALSE, sep='\t')
dat2 = read.delim("ps_Del20_challenge.txt", header = FALSE, sep='\t')
head(dat1)
# combine two ps together
dat_deltap <- data.frame(Contrast = factor(rep(c("Del19","Del20"), each=length(dat1$V5))), Delta_p = c(dat1$V5, dat2$V5))
# start plot
cdat <- ddply(dat_deltap, "Contrast", summarise, rating.mean=mean(Delta_p))
cdat
tiff("delta_p.tiff", units="in", width=8, height=5, res=300)
# Density plots with semi-transparent fill (only focus on count https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1)
ggplot(dat_deltap, aes(x=Delta_p, fill=Contrast)) + 
  geom_histogram(alpha=.6, bins = 80, position = "identity") +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=Contrast), linetype="dashed", size=1, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  xlab(expression(Delta~italic(p))) + ylab("Frequency")+
  guides(fill=guide_legend(title="Contrast"))+
  theme_bw() +
  theme(text = element_text(size=16,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
dev.off()

# load the ps datasets
dat1 = read.delim('ps_Del19_challenge.txt', header = FALSE, sep='\t')
dat2 = read.delim("ps_Del20_challenge.txt", header = FALSE, sep='\t')
head(dat1)
length(dat1$V6)
# combine two ps together
dat_ps <- data.frame(Contrast = factor(rep(c("Del19","Del20"), each=2032113)), ps = c(dat1$V6, dat2$V6))
# start plot
cdat <- ddply(dat_ps, "Contrast", summarise, rating.mean=mean(ps))
cdat
tiff("ps.tiff", units="in", width=8, height=5, res=300)
# Density plots with semi-transparent fill (only focus on count https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1)
ggplot(dat_ps, aes(x=ps, fill=Contrast)) + 
  geom_histogram(alpha=.6, bins=50, binwidth=0.055, position = "identity") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  xlab(expression(italic(p)~values)) + ylab("Frequency")+
  guides(fill=guide_legend(title="Contrast"))+
  theme_bw() +
  theme(text = element_text(size=16,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
dev.off()

# A test data from global angsd calling (1x but not shared SNPs)
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_05_saf_maf_by_pop/test")
CHR19 = read.delim('CHR19_all_minq20_minmq30_CV30_masked.mafs', header = TRUE, sep='\t')
REF19 = read.delim('REF19_all_minq20_minmq30_CV30_masked.mafs', header = TRUE, sep='\t')
CHR20 = read.delim('CHR20_all_minq20_minmq30_CV30_masked.mafs', header = TRUE, sep='\t')
REF20 = read.delim('REF20_all_minq20_minmq30_CV30_masked.mafs', header = TRUE, sep='\t')
delta1 = CHR19$knownEM-REF19$knownEM
delta2 = CHR20$knownEM-REF20$knownEM
dat_deltap <- data.frame(Contrast = factor(rep(c("Del19","Del20"), each=length(delta1))), Delta_p = c(delta1, delta2))
# start plot
cdat <- ddply(dat_deltap, "Contrast", summarise, rating.mean=mean(Delta_p))
cdat
tiff("delta_p.tiff", units="in", width=8, height=5, res=300)
# Density plots with semi-transparent fill (only focus on count https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1)
ggplot(dat_deltap, aes(x=Delta_p, fill=Contrast)) + 
  geom_histogram(alpha=.6, bins = 80, position = "identity") +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=Contrast), linetype="dashed", size=1, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  xlab(expression(Delta~italic(p))) + ylab("Frequency")+
  guides(fill=guide_legend(title="Contrast"))+
  theme_bw() +
  theme(text = element_text(size=16,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
dev.off()


# A test data 
set.seed(1234)
testA = rchisq(1000, df = 1)
testA = (testA - min(testA)) / (max(testA) - min(testA))
testA.adj = p.adjust(testA, method = 'BH')
testA.adj[testA.adj<0.1]

testB = rchisq(1000, df = 4)
testB = (testB - min(testB)) / (max(testB) - min(testB))
testB.adj = p.adjust(testB, method = 'BH')
testB.adj[testB.adj<0.1]

test <- data.frame(cond = factor(rep(c("A","B"), each=1000)), 
                  ps = c(testA, testB))

# Density plots with semi-transparent fill
ggplot(data = test) +
  geom_histogram(alpha=.6, mapping = aes(x = ps, fill = cond), position = "identity", bins=80)




################ for DelBay19_dataset ################
setwd("~/Documents/Ryan_workplace/DelBay19_adult/11_SGS")
# load the ps datasets
dat1 = read.delim('ps_Del19_challenge.txt', header = FALSE, sep='\t')
dat2 = read.delim("ps_Del19_HC_NB.txt", header = FALSE, sep='\t')
head(dat1)
# combine two ps together
dat_deltap <- data.frame(Contrast = factor(rep(c("Del19","HC_NB"), each=length(dat1$V5))), Delta_p = c(dat1$V5, dat2$V5))
# start plot
cdat <- ddply(dat_deltap, "Contrast", summarise, rating.mean=mean(Delta_p))
cdat
tiff("delta_p.tiff", units="in", width=8, height=5, res=300)
# Density plots with semi-transparent fill (only focus on count https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1)
ggplot(dat_deltap, aes(x=Delta_p, fill=Contrast)) + 
  geom_histogram(alpha=.6, bins=50, binwidth=0.055, position = "identity") +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=Contrast), linetype="dashed", size=1, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  xlab(expression(Delta~italic(p))) + ylab("Frequency")+
  guides(fill=guide_legend(title="Contrast"))+
  theme_bw() +
  theme(text = element_text(size=16,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
dev.off()

# load the ps datasets
dat1 = read.delim('ps_Del19_challenge.txt', header = FALSE, sep='\t')
dat2 = read.delim("ps_Del19_HC_NB.txt", header = FALSE, sep='\t')
dat3 = read.delim("ps_Del19_HC_SR.txt", header = FALSE, sep='\t')
head(dat1)
length(dat1$V3)
# combine two ps together
dat_ps <- data.frame(Contrast = factor(rep(c("Del19","HC_NB", "HC_SR"), each=length(dat1$V3))), ps = c(dat1$V3, dat2$V3, dat3$V3))
# start plot
cdat <- ddply(dat_ps, "Contrast", summarise, rating.mean=mean(ps))
cdat
tiff("ps.tiff", units="in", width=8, height=5, res=300)
# Density plots with semi-transparent fill (only focus on count https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1)
ggplot(dat_ps, aes(x=ps, fill=Contrast)) + 
  geom_histogram(alpha=.3, bins=50, binwidth=0.055, position = "identity") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  xlab(expression(italic(p)~values)) + ylab("Frequency")+
  guides(fill=guide_legend(title="Contrast"))+
  theme_bw() +
  theme(text = element_text(size=16,family="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
dev.off()
