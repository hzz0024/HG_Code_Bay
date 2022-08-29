require(graphics)
require(dgof)

x <- runif(50)
y <- runif(30)

ks.test(x, y)
barplot(y[-1])
barplot(x[-1])

rm(list=ls())
#############################
# format the Fisher outlier #
#############################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/25_KS_sfs")
format_outlier <- function(headname){
  #headname = "18_SGS_HC_NB_FDR_outlier.list"
  dat <- read.table(headname, header = T)
  # do formatting for each chromosome
  #chr = c('NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')
  #for (j in 1:10){ 
  #  dat$V1[dat$V1 == j] = chr[j]
  #}
  dat$snp <- paste0(dat$chr,'_',dat$pos)
  ref <- read.table("All_sfs_all_sites_noparalogs.snplist.txt", header = F)
  outlier_list <- ref[which(paste0(ref$V1,'_',ref$V2) %in% dat$snp), ]
  print(paste0("number of outliers:", dim(outlier_list[1])))
  neutral_list <- ref[-which(paste0(ref$V1,'_',ref$V2) %in% dat$snp), ]
  print(paste0("number of neutral snps:", dim(neutral_list[1])))
  write.table(outlier_list, paste0(strsplit(headname, split = ".list")[[1]][1], ".outlier.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  write.table(neutral_list, paste0(strsplit(headname, split = ".list")[[1]][1], ".neutral.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_outlier("18_SGS_HC_NB_FDR_outlier.list")
format_outlier("19_SGS_HC_NB_FDR_outlier.list")
format_outlier("21_SGS_HC_NB_FDR_outlier.list")

#############################
# plot an example           #
#############################

par(mar=c(2,4,1,4))
par(mfrow=c(5,2))

plot_sfs1 <- function(head){
  sfs<-scan(paste0(head, ".sfs"))
  #function to normalize
  norm <- function(x) x/sum(x)
  #the variable categories of the sfs
  sfs<-norm(sfs[-c(1,length(sfs))]) 
  barplot(sfs, xlab=NULL, main=paste0(head),
          names=1:length(sfs),ylab="Proportions",col='blue')
}
plot_sfs2 <- function(head){
  sfs<-scan(paste0(head, ".sfs"))
  #function to normalize
  norm <- function(x) x/sum(x)
  #the variable categories of the sfs
  sfs<-norm(sfs[-c(1,length(sfs))]) 
  barplot(sfs, xlab=NULL, main=paste0(head),
          names=1:length(sfs),ylab="Proportions",col='orange')
}
plot_sfs3 <- function(head){
  sfs<-scan(paste0(head, ".sfs"))
  #function to normalize
  norm <- function(x) x/sum(x)
  #the variable categories of the sfs
  sfs<-norm(sfs[-c(1,length(sfs))]) 
  barplot(sfs, xlab="Number of minor alleles", main=paste0(head),
          names=1:length(sfs),ylab="Proportions",col='blue')
}
plot_sfs4 <- function(head){
  sfs<-scan(paste0(head, ".sfs"))
  #function to normalize
  norm <- function(x) x/sum(x)
  #the variable categories of the sfs
  sfs<-norm(sfs[-c(1,length(sfs))]) 
  barplot(sfs, xlab="Number of minor alleles", main=paste0(head),
          names=1:length(sfs),ylab="Proportions",col='orange')
}

plot_all <-function(year){
  plot_sfs1(paste0("HC_", year, ".neutral"))
  plot_sfs2(paste0("HC_", year, ".outlier"))
  plot_sfs1(paste0("ARN_", year, ".neutral"))
  plot_sfs2(paste0("ARN_", year, ".outlier"))
  plot_sfs1(paste0("COH_", year, ".neutral"))
  plot_sfs2(paste0("COH_", year, ".outlier"))
  plot_sfs1(paste0("SR_", year, ".neutral"))
  plot_sfs2(paste0("SR_", year, ".outlier"))
  plot_sfs3(paste0("NB_", year, ".neutral"))
  plot_sfs4(paste0("NB_", year, ".outlier"))
  #graph2ppt(file=past0("WILD_",year),width=10,height=6)
}

par(mar=c(4,4,4,2))
par(mfrow=c(5,2))
plot_all(18)
plot_all(19)
plot_all(21)
#export as 1400 1050

#############################
# Extract the  outlier      #
#############################

extract_snp <- function(outlier, full_list){
  #outlier = "challenge_FDR_outlier.list"
  dat <- read.table(outlier, header = F)
  # do formatting for each chromosome
  chr = c('NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')
  for (j in 1:10){ 
    dat$V1[dat$V1 == j] = chr[j]
  }
  dat$id <- paste0(dat$V1,'_',dat$V2)

  #full_list = "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt" this is initially used to build the sfs with realSFS. Now I only focused on the maf values
  # full <- read.table(full_list, header = F)
  # full$id <- paste0(full$V1,'_',full$V2)
  # angsd_snplist <- full[which(full$id %in% dat$id), ][1:4]
  # write.table(angsd_snplist, paste0(strsplit(outlier, split = ".list")[[1]][1], ".snplist.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  
  #full_list = "REF19_all_minq20_minmq30_CV30_masked.mafs"
  full <- read.table(full_list, header = T)
  full$id <- paste0(full$chromo,'_',full$position)
  angsd_snplist <- full[which(full$id %in% dat$id), ]$knownEM
  write.table(angsd_snplist, paste0(strsplit(outlier, split = ".list")[[1]][1], ".mafs.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  return(angsd_snplist)
}

outlier <- extract_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "REF19_all_minq20_minmq30_CV30_masked.mafs")
outlier <- extract_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "NB_all_minq20_minmq30_CV30_masked.mafs")
outlier <- extract_snp("challenge_FDR_outlier.list", "REF19_all_minq20_minmq30_CV30_masked.mafs")
outlier <- extract_snp("HC_NB_FDR_outlier.list", "NB_all_minq20_minmq30_CV30_masked.mafs")
# extract_snp("REF19_CHR19_NB_HC_out_0.05_fish.txt", "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
# extract_snp("challenge_FDR_outlier.list", "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
# extract_snp("HC_NB_FDR_outlier.list", "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")

#############################
# Select random SNPs        #
#############################

random_snp <- function(outlier, full_list){
  #outlier = "challenge_FDR_outlier.list"
  dat <- read.table(outlier, header = F)
  # do formatting for each chromosome
  chr = c('NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')
  for (j in 1:10){ 
    dat$V1[dat$V1 == j] = chr[j]
  }
  dat$id <- paste0(dat$V1,'_',dat$V2)
  
  #full_list = "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt" this is initially used to build the sfs with realSFS. Now I only focused on the maf values
  # full <- read.table(full_list, header = F)
  # full$id <- paste0(full$V1,'_',full$V2)
  # angsd_snplist <- full[which(full$id %in% dat$id), ][1:4]
  # write.table(angsd_snplist, paste0(strsplit(outlier, split = ".list")[[1]][1], ".snplist.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  
  #full_list = "REF19_all_minq20_minmq30_CV30_masked.mafs"
  full <- read.table(full_list, header = T)
  full$id <- paste0(full$chromo,'_',full$position)
  angsd_snplist <- full[sample(dim(full)[1], dim(dat)[1]), ]$knownEM
  write.table(angsd_snplist, paste0(strsplit(outlier, split = ".list")[[1]][1], ".random.mafs.txt"), row.names=F, col.names = F, quote=F, sep="\t")
  return(angsd_snplist)
}

random <- random_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "REF19_all_minq20_minmq30_CV30_masked.mafs")
random <- random_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "NB_all_minq20_minmq30_CV30_masked.mafs")

random <- random_snp("challenge_FDR_outlier.list", "REF19_all_minq20_minmq30_CV30_masked.mafs")
random <- random_snp("HC_NB_FDR_outlier.list", "NB_all_minq20_minmq30_CV30_masked.mafs")
# random_snp("challenge_FDR_outlier.list", "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
# random_snp("HC_NB_FDR_outlier.list", "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")

#############################
# p-value                   #
#############################
library(plyr)
df_total = as.data.frame(c(outlier, random))
df_total$group = rep(c("outlier", "random"), c(length(outlier), length(random)))
colnames(df_total)=c('maf', 'group')

mu <- ddply(df_total, "group", summarise, grp.mean=mean(maf))
head(mu)

outlier <- extract_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "REF19_all_minq20_minmq30_CV30_masked.mafs")

ps = c()
for(i in seq(100)){
  tmp = mean(random_snp("REF19_CHR19_NB_HC_out_0.05_fish.snplist.txt", "NB_all_minq20_minmq30_CV30_masked.mafs"))
  ps = c(ps, tmp)
}

p = (sum(ps > mean(outlier)))/100

#############################
# Perform ks.test           #
#############################
require(graphics)
require(dgof)
source("signed-ks-test.R")

ks.test(outlier, random)
barplot(outlier)
barplot(random)

#############################
# hisgram plot              #
#############################
library(plyr)
df_total = as.data.frame(c(outlier, random))
df_total$group = rep(c("outlier", "random"), c(length(outlier), length(random)))
colnames(df_total)=c('maf', 'group')

mu <- ddply(df_total, "group", summarise, grp.mean=mean(maf))
head(mu)

ggplot(df_total, aes(x=maf, color=group, fill=group)) +
  geom_histogram(alpha=0.5, bins = 30, fill= "#FFFAFA", position="identity")+
  #geom_vline(aes(xintercept=mean(random_mean)), color="grey",linetype="dashed")+
  labs(x="Allele frequency", y = "Count")+
  theme_classic() +
  theme(text = element_text(size = 15)) +
  theme(legend.position="right")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")

#############################
# plot an example           #
#############################
sfs<-scan("COH_18.neutral.sfs")
barplot(sfs[-1])
barplot(sfs[-c(1,length(sfs))])

par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
#function to normalize
norm <- function(x) x/sum(x)
#read data
sfs <- (scan("challenge_FDR_outlier.sfs"))
sfs <- (scan("challenge_FDR_random.sfs"))
#the variable categories of the sfs
sfs<-norm(sfs[-c(1,length(sfs))]) 

barplot(sfs, xlab="Number of minor alleles",
        names=1:length(sfs),ylab="Proportions",col='blue')

