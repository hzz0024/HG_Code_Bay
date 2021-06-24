setwd("~/Documents/Ryan_workplace/DelBay_adult/14_PRF")

library(tidyr)
library(dplyr)
library(ggplot2)

df <- read.csv("Del19_p001_summary.csv")
df$ntree <- factor(df$ntree)
head(df)
ggplot(df, aes(x=Qutantile, y=accuracy)) + 
  geom_line(aes(color = ntree)) +
  xlab("Percentile of top SNPs") + ylab("Accuracy")

head(df)
ggplot(df, aes(x=Qutantile, y=accuracy)) + 
  geom_line(aes(color = ntree)) + xlim(0, 0.1)+
  xlab("Percentile of top SNPs") + ylab("Accuracy")

df <- data.frame(x = c(0, 5, 10, 15), y = c(2.2, 3.8, 4.6, 7.6), z = c(4.5, 6.8, 9.3, 10.5))
ggplot(df, aes(x)) + 
  geom_line(aes(y = y, colour = "y")) + 
  geom_line(aes(y = z, colour = "z")) +


setwd("~/Documents/Ryan_workplace/DelBay_adult/14_PRF")
df <- read.csv("Del20_FDR005_2x_summary.csv")

df$ntree <- factor(df$ntree)
head(df)
ggplot(df, aes(x=SNPs, y=accuracy)) + 
  geom_line(aes(color = ntree)) +
  xlab("Top SNP numbers") + ylab("Accuracy") +
  geom_vline(xintercept=182, linetype="dashed", color = "red", size=0.5)

############################################
#### Delta-p - importance recorrelation ####
############################################
setwd("/Volumes/cornell/beagle_for_PRF")
name = "ps_Del20_challenge.txt" # 2-side p-value
dat1 = read.delim(name, header = FALSE, sep='\t')
dat1$id = paste0(dat1$V1,'_',dat1$V2)
dat1[with(dat1, order(id)),]
head(dat1)

name = "toy_test_importance.txt" # 2-side p-value
dat2 = read.delim(name, header = FALSE, sep='\t')
dat = dat2[with(dat2, order(V1)),]
head(dat)
# delta_p 
delta_p <- dat1[which(dat1$id %in% dat2$V1),][,1:5]
delta_p$id = dat1[which(dat1$id %in% dat2$V1),]$id
delta_p = delta_p[with(delta_p, order(id)),]
# importance

dat$V1 == delta_p$id

library("ggpubr")
cor(dat$V2, abs(delta_p$V5), method = 'spearman')
qqplot(abs(delta_p$V5), dat$V2*25, main='364')
abline(c(0,1))

name = "Del20_FDR_outlier.list" # 2-side p-value
dat0 = read.delim(name, header = TRUE, sep='\t')
dat0$id = paste0(dat0$chr,'_',dat0$pos)
dat0[with(dat0, order(id)),]

color = delta_p$id %in% dat0$id
#color[color==FALSE] = 2
#color[color==TRUE] = 0
color[color==FALSE] = 2
color[color==TRUE] = 1

# Fst distribution for 182 Del20 SGS outliers
hist(abs(delta_p$V5)[delta_p$id %in% dat0$id])
# qqplot for delta_p vs importance 
#qqplot(abs(delta_p$V5), dat$V2*25, main='364',col=color)
# plot delta_p (absolute value) vs importance 
plot(abs(delta_p$V5), dat$V2, main='364',col=color)
abline(lm(dat$V2~abs(delta_p$V5)), col="blue") # regression line (y~x)
# plot delta_p (actual values) vs importance
plot((delta_p$V5), dat$V2, main='364',col=color)
#abline(lm(dat$V2~(delta_p$V5)), col="blue") # regression line (y~x)


setwd("/Volumes/cornell/beagle_for_PRF")
name = "Del20_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.mafs" # 2-side p-value
dat1 = read.delim(name, header = TRUE, sep='\t')
dat1$id = paste0(dat1$chromo,'_',dat1$position)
dat1[with(dat1, order(id)),]
head(dat1)

name = "toy_test_importance.txt" # 2-side p-value
dat2 = read.delim(name, header = FALSE, sep='\t')
dat = dat2[with(dat2, order(V1)),]
head(dat)
# delta_p 
delta_p <- dat1[which(dat1$id %in% dat2$V1),]
delta_p$id = dat1[which(dat1$id %in% dat2$V1),]$id
delta_p = delta_p[with(delta_p, order(id)),]
# importance

dat$V1 == delta_p$id

library("ggpubr")
cor(dat$V2, abs(delta_p$V5), method = 'spearman')
qqplot(abs(delta_p$V5), dat$V2*25, main='364')
abline(c(0,1))

name = "Del20_FDR_outlier.list" # 2-side p-value
dat0 = read.delim(name, header = TRUE, sep='\t')
dat0$id = paste0(dat0$chr,'_',dat0$pos)
dat0[with(dat0, order(id)),]

color = delta_p$id %in% dat0$id
#color[color==FALSE] = 2
#color[color==TRUE] = 0
color[color==FALSE] = 2
color[color==TRUE] = 1

# Fst distribution for 182 Del20 SGS outliers
hist(abs(delta_p$knownEM)[delta_p$id %in% dat$V1])
# qqplot for delta_p vs importance 
#qqplot(abs(delta_p$V5), dat$V2*25, main='364',col=color)
# plot delta_p (absolute value) vs importance 
plot(dat_plot$maf, dat_plot$imp, main='364',col=color)
abline(lm(dat$V2[dat$V1 %in% delta_p$id]~abs(delta_p$knownEM)), col="blue") # regression line (y~x)
# plot delta_p (actual values) vs importance
plot((delta_p$knownEM), dat$V2[dat$V1 %in% delta_p$id], main='364',col=color)
#abline(lm(dat$V2~(delta_p$V5)), col="blue") # regression line (y~x)

dat$V1[dat$V1 %in% delta_p$id]
dat_plot <- data.frame(id1 = dat$V1[dat$V1 %in% delta_p$id], imp = dat$V2[dat$V1 %in% delta_p$id],  id2 = delta_p$id[delta_p$id %in% dat$V1], maf =delta_p$knownEM[delta_p$id %in% dat$V1])



tmp=cbind(dat$V2, delta_p$V5)
tmp[with(dat, order(V2)),]



name = "Del20_FDR_outlier.list" 
dat1 = read.delim(name, header = TRUE, sep='\t')
dat1$id = paste0(dat1$chr,'_',dat1$pos)
dat1[with(dat1, order(id)),]

name = "toy_test_importance.txt" # 2-side p-value
dat2 = read.delim(name, header = FALSE, sep='\t')
#dat = dat2[with(dat2, order(V1)),]

imp = dat2[which(dat2$V1 %in% dat1$id),][,1:2]
imp = imp[with(imp, order(V1)),]

imp$V1 == dat1$id
cor(imp$V2, abs(dat1$deltap), method = 'spearman')
