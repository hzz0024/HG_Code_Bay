#############################################
## reorder the ps file from SGS test       ## 
############################################# 

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/no_shared/indep_SGS_ps")
########### be really careful about the issue from order. clear everthing in the data before doing that.
pname = "ps_CHR19_REF19.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
dat = dat[order(dat[,1], dat[,2]),]
write.table(dat, "./ps_CHR19_REF19_reorder.txt", row.names=F, col.names=F, quote=F, sep="\t")

####################################
##########  output snp list ########
####################################
# output global snp list
df = data.frame(dat$V1, dat$V2)
write.table(df, "./SNP_list/SGS_global_snp_list.txt", row.names=F, col.names = F, quote=F, sep="\t")
# output outliers with FDR < 0.05
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'delta_p', 'ps', 'raw_candidates', 'adj')
dat1 <- dat[which(dat$adj< 0.05),]
dat1$id <- paste0(dat1$chromo,'_',dat1$position)
dim(dat1)
dat1 = dat1[with(dat1, order(chromo, position)),]
df1 <- data.frame(dat1$chromo, dat1$position)
write.table(df1, "./SNP_list/CHR19_REF19_SGS_FDR05.txt", row.names=F, col.names = F, quote=F, sep="\t")

#############################################
## Check angular deltap in combined p SGS ## 
############################################# 

check_direct <- function(pname){
  #pname = "./REF19_CHR19_REF20_CHR20_out_0.05_fish.txt"
  dat = read.delim(pname, header = TRUE, sep='\t')
  dp1 = dat$p_chr -	dat$p_ref
  dp2 = dat$p_HC -	dat$p_NB
  share_all = sum(sign(as.numeric(dp1)) == sign(as.numeric(dp2)))
  message("number of SNPs with same directionality is ", share_all, " out of ", dim(dat)[1], " (", share_all/dim(dat)[1], ")")
}

check_direct("REF19_CHR19_REF20_CHR20_out_0.05_fish.txt")
check_direct("REF19_CHR19_REF20_CHR20_out_0.01_fish.txt")
check_direct("REF19_CHR19_REF20_CHR20_out_all_fish.txt")
check_direct("SR_REF19_ARN_COH_out_0.05_fish.txt")
check_direct("SR_REF19_ARN_COH_out_0.01_fish.txt")
check_direct("SR_REF19_ARN_COH_out_all_fish.txt")
check_direct("REF19_CHR19_HC_NB_out_0.01_fish.txt")
check_direct("REF19_CHR19_HC_NB_out_0.05_fish.txt")
check_direct("REF19_CHR19_HC_NB_out_all_fish.txt")

#############################################
## Check angular deltap in independent SGS ## 
############################################# 

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/indep_SGS_angular_ps")

# test
# load reference file with header in it
ref = 'REF19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')

# load challenge file with header in it
chr = 'CHR19_minmapq30_minq20_CV30_masked_noinvers_shared_sites.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')

AFS <- function(theta,N)
{
  Pp = NULL
  for (j in 1:(2*N-1))
  {Pp = c(Pp,theta*((1/j)+(1/(2*N-j))))}
  Pp = Pp*(1/sum(Pp))
}

Dfreq <- function(X,N,theta){
  theta = as.numeric(theta)
  # p is the vector that stores all the frequencies that the
  # allele can take in a sample of size N
  p=(1:(2*N-1))/(2*N)
  # we get P(X|p) over all possible values of p within the pool
  # this is equivalent to exp(lfactorial(N)-(lfactorial(N - X)+lfactorial(X))) * p^X * (1-p)^(N-K)
  PXp = dbinom(X,2*N,p)
  # we get the probability density of allele counts P(p):
  Pp = AFS(theta,N)
  # and P(X) over [1,2N-1] is simply P(X) = 1/(2N-1)
  PX = 1/(2*N-1)
  # So for a given value of X observed in a sample of size N
  # we get P(p|X) the probability density of p in the gene pool
  PpX = NULL
  for(k in 1:(2*N-1))
  {PpX=c(PpX,(PXp[k]*Pp[k])/PX)}
  # and we scale it to one
  PpX = PpX*(1/sum(PpX))
}

Ddelta <- function(X, N, K, theta, nreps)
{
  poolfreq = colSums((1:(2*N-1))*rmultinom(nreps,1:(2*N-1),prob=matrix(Dfreq(X,N,theta),1)))/(2*N)
  samplefreq = rbinom(nreps,2*K,poolfreq)/(2*K)
  #delta = samplefreq - (X/(2*N))
  delta = abs(2*asin(sqrt(samplefreq)) - 2*asin(sqrt((X/(2*N)))))/pi
}

nreps = 100000
i=1
ref_n = dat_ref$nInd[i]
ref_k = round(dat_ref$knownEM[i]*dat_ref$nInd[i]*2)
ch_n = dat_ch$nInd[i]
theta = 0.015965112 # replace the local theta with global theta (i.e. 2019 reference theta)
obs_delta=abs(2*asin(sqrt(dat_ch$knownEM[i])) - 2*asin(sqrt(dat_ref$knownEM[i])))/pi
delta_ps = Ddelta(X=ref_k, N=ref_n, K=ch_n, theta=theta, nreps=nreps)
p_value <- (length(delta_ps[abs(delta_ps) >= abs (obs_delta)]))/(length(delta_ps))
# plot
hist(delta_ps, main="Example 1", 
     xlab="Angular transformed delta_p" )
abline(v = obs_delta, col="red", lwd=3, lty=2)

nreps = 100000
i=292
ref_n = dat_ref$nInd[i]
ref_k = round(dat_ref$knownEM[i]*dat_ref$nInd[i]*2)
ch_n = dat_ch$nInd[i]
theta = 0.015965112 # replace the local theta with global theta (i.e. 2019 reference theta)
obs_delta=abs(2*asin(sqrt(dat_ch$knownEM[i])) - 2*asin(sqrt(dat_ref$knownEM[i])))/pi
delta_ps = Ddelta(X=ref_k, N=ref_n, K=ch_n, theta=theta, nreps=nreps)
p_value <- (length(delta_ps[abs(delta_ps) >= abs (obs_delta)]))/(length(delta_ps))
# plot
hist(delta_ps, main="Example 2", 
     xlab="Angular transformed delta_p" )
abline(v = obs_delta, col="red", lwd=3, lty=2)

#############################################
# A simple comp between angular vs original #
#############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/indep_SGS_angular_vs_orginals")
pname = "original_ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')

cnt_outlier <- length(dat$adj[which(dat$adj< 0.05)])
angular_outlier <- dat[which(dat$adj< 0.05),]
id1 <- paste0(angular_outlier$chromo,'_',angular_outlier$position)

neg <- length(angular_outlier$adj[which(angular_outlier$delta_p<0)])
pos <-  length(angular_outlier$adj[which(angular_outlier$delta_p>=0)])
message("number of outliers with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

neg <- length(dat$delta_p[which(dat$delta_p<0)])
pos <- length(dat$delta_p[which(dat$delta_p>=0)])
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

dat$adj = -log(dat$adj)

source("manhattan.R")
jpeg("Mahattan_adj_comp.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(2,1))
manhattan(chr="chromo",bp="position",p="adj", snp = "SNP", dat, highlight1 = id1 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 5),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log(FDR) using original delta_P", cex.lab=1.5) 

pname = "angular_ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')

cnt_outlier <- length(dat$adj[which(dat$adj< 0.05)])
angular_outlier <- dat[which(dat$adj< 0.05),]
id2 <- paste0(angular_outlier$chromo,'_',angular_outlier$position)

neg <- length(angular_outlier$adj[which(angular_outlier$delta_p<0)])
pos <-  length(angular_outlier$adj[which(angular_outlier$delta_p>=0)])
message("number of outliers with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

neg <- length(dat$delta_p[which(dat$delta_p<0)])
pos <- length(dat$delta_p[which(dat$delta_p>=0)])
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

dat$adj = -log(dat$adj)

manhattan(chr="chromo",bp="position",p="adj", snp = "SNP", dat, highlight1 = id2 , logp=FALSE, cex.axis = 1.2, ylim = c(0, 5),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="-log(FDR) using Angular delta_P", cex.lab=1.5) 

length(intersect(id1,id2))
# 1605
# > length(id2) # Angular
# [1] 3550
# > length(id1) # Original
# [1] 2073

dev.off()

############################################
## Check SNP direction in independent SGS ## 
############################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/indep_SGS_angular_ps")
# check the delta_p directionality
check_dir <- function(pname){
  #pname = "ps_Del19_challenge.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  message("number of SNPs with positive D' is ", length(dat$V5[which(dat$V5>=0)]))
  dat$delta_p <- dat$V3-dat$V4
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
  # process outliers
  cnt_outlier <- length(dat$adj[which(dat$adj< 0.05)])
  angular_outlier <- dat[which(dat$adj< 0.05),]
  angular_outlier$id <- paste0(angular_outlier$chromo,'_',angular_outlier$position)
  # plot
  # jpeg("Angular_transformation_before_after.jpg", width = 12, height = 8, units = 'in', res = 300)
  # par(mfrow=c(2,1))
  # # plot the af distribution in 2019 Ref
  # hist(angular_outlier$p0, main="Allele frequency in 2019 Ref. before transformation", 
  #      xlab="Allele frequency distribution in 2019 Ref." )
  # 
  # hist(2*asin(sqrt(angular_outlier$p0))/pi, main="Allele frequency in 2019 Ref. after Angular transformation", 
  #      xlab="Angular allele frequency distribution in 2019 Ref." )
  # dev.off()
  neg <- length(angular_outlier$adj[which(angular_outlier$delta_p<0)])
  pos <-  length(angular_outlier$adj[which(angular_outlier$delta_p>=0)])
  message("number of angular outliers with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")
  
  neg <- length(dat$delta_p[which(dat$delta_p<0)])
  pos <- length(dat$delta_p[which(dat$delta_p>=0)])
  message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")
  
  p<-ggplot(data=dat1, aes(ps)) + 
    geom_density(alpha = 0.2, size = 1.2) + theme(legend.position = "right") +
    #geom_histogram(bins = 40, color="black", fill="white")+
    scale_x_continuous(expression(italic(p)~"-"~values)) + 
    ggtitle(paste0("Genome-wide p-values in ", pname))+
    theme(axis.text=element_text(size=14),
          text = element_text(size=14), # add family="times" in macOS for time font
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "white", size=0.5),
          plot.margin=unit(c(1,1,1.5,1.2),"cm")) 
  print(p)
  #hist(angular_outlier$adj, main=paste0("p-values distribution in ", pname), 
  #     xlab=expression(italic(p)~"-"~values))
}

check_dir("ps_Del19_challenge.txt")
check_dir("ps_Del20_challenge.txt")
check_dir("ps_Del19_HC_NB.txt")
check_dir("ps_SR_REF19.txt")
check_dir("ps_ARN_COH.txt")

##############################################
## plot pvalue distribution among contrasts ## 
##############################################

pname1 = "ps_Del19_challenge.txt"
dat1 = read.delim(pname1, header = FALSE, sep='\t')
message("number of SNPs with positive D' is ", length(dat1$V5[which(dat1$V5>=0)]))
dat1$delta_p <- dat1$V3-dat1$V4
dat1$adj = p.adjust(dat1$V6, method = 'BH')
colnames(dat1)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
dat1 <- dat1[which(dat1$adj< 0.05),]
dat1$id <- paste0(dat1$chromo,'_',dat1$position)
df1 <- data.frame(dat1$id, dat1$ps, "2019 Challenge")
colnames(df1)=c('id', 'ps', 'group')

pname2 = "ps_Del19_HC_NB.txt"
dat2 = read.delim(pname2, header = FALSE, sep='\t')
dat2$id <- paste0(dat2$V1,'_',dat2$V2)
dat2 = dat2[which(dat2$id %in% dat1$id),]
df2 <- data.frame(dat2$id, dat2$V6, "Wild HC vs NB")
colnames(df2)=c('id', 'ps', 'group')

#pname3 = "ps_SR_REF19.txt"
pname3 = "ps_Del19_SR_NB.txt"
dat3 = read.delim(pname3, header = FALSE, sep='\t')
dat3$id <- paste0(dat3$V1,'_',dat3$V2)
dat3 = dat3[which(dat3$id %in% dat1$id),]
#df3 <- data.frame(dat3$id, dat3$V6, "Wild SR vs 2019 Ref")
df3 <- data.frame(dat3$id, dat3$V6, "Wild SR vs NB")
colnames(df3)=c('id', 'ps', 'group')

pname4 = "ps_Del20_challenge.txt"
dat4 = read.delim(pname4, header = FALSE, sep='\t')
dat4$id <- paste0(dat4$V1,'_',dat4$V2)
dat4 = dat4[which(dat4$id %in% dat1$id),]
df4 <- data.frame(dat4$id, dat4$V6, "2020 Challenge")
colnames(df4)=c('id', 'ps', 'group')

#pname5 = "ps_ARN_COH.txt"
pname5 = "ps_Del19_HC_SR.txt"
dat5 = read.delim(pname5, header = FALSE, sep='\t')
dat5$id <- paste0(dat5$V1,'_',dat5$V2)
dat5 = dat5[which(dat5$id %in% dat1$id),]
#df5 <- data.frame(dat5$id, dat5$V6, "Wild ARN vs COH")
df5 <- data.frame(dat5$id, dat5$V6, "Wild HC vs SR")
colnames(df5)=c('id', 'ps', 'group')

df <- rbind(df1, df2, df3, df4, df5)
jpeg("Indep_SGS_p-value_comp.jpg", width = 12, height = 8, units = 'in', res = 300)
p<-ggplot(data=df, aes(x=ps, color=group)) + 
  #geom_density(alpha = 0.2, size = 1.2) + theme(legend.position = "right") +
  geom_freqpoly(bins = 20)+
  scale_x_continuous(expression(italic(p)~"-"~values)) + 
  scale_y_continuous("Count") +
  ggtitle(paste0("P-values (before FDR correction) distribution for outliers identified from \n2019 Challenge, and corrsponding p-values in other contrasts"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "white", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 
p
dev.off()










par(mfrow=c(2,1))
jpeg("SNP_FDR05_Del19.jpg", width = 12, height = 8, units = 'in', res = 300)
p<-ggplot(data=df1, aes(ps)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(SNP~candidates~Delta~italic(p)~"in"~'2019'~challenge)) + 
  ggtitle("Candidates with positive changes is \n1805 out of 2185 (0.83)")+
  #ggtitle("Candidates with positive changes is \n441 out of 540 (0.82)")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname1 = "ps_Del19_HC_NB.txt"
dat1 = read.delim(pname1, header = FALSE, sep='\t')
dat1$id1 = paste0(dat1$V1,'_',dat1$V2)
dat1_outlier = dat1[which(dat1$id1 %in% outlier_id),]
neg <- length(dat1_outlier$V5[which(dat1_outlier$V5<0)])
pos <- length(dat1_outlier$V5[which(dat1_outlier$V5>0)])
# 1042/(2185)
# [1] 0.4768
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p1<-ggplot(data=dat1_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~wild~transect~"("~HC~-~NB~")")) + 
  ggtitle("Same SNPs in HC - NB but ratio \nis 1042 out of 2185 (0.48) ")+
  #ggtitle("Same SNPs in HC - NB but ratio \nis 248 out of 540 (0.46) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname2 = "ps_Del19_SR_NB.txt"
dat2 = read.delim(pname2, header = FALSE, sep='\t')
dat2$id1 = paste0(dat2$V1,'_',dat2$V2)
dat2_outlier = dat2[which(dat2$id1 %in% outlier_id),]
neg <- length(dat2_outlier$V5[which(dat2_outlier$V5<0)])
pos <- length(dat2_outlier$V5[which(dat2_outlier$V5>0)])
# 1109/(2185)
# [1] 0.50.8
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p2<-ggplot(data=dat2_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~wild~transect~"("~SR~-~NB~")")) + 
  ggtitle("Same SNPs in SR - NB but ratio \nis 1057 out of 2185 (0.48) ")+
  #ggtitle("Same SNPs in SR - NB but ratio \nis 299 out of 540 (0.55) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

# pname2 = "ps_Del19_ARN_COH.txt"
# dat2 = read.delim(pname2, header = FALSE, sep='\t')
# dat2$id1 = paste0(dat2$V1,'_',dat2$V2)
# dat2_outlier = dat2[which(dat2$id1 %in% outlier_id),]
# neg <- length(dat2_outlier$V5[which(dat2_outlier$V5<0)])
# pos <- length(dat2_outlier$V5[which(dat2_outlier$V5>0)])
# # 1109/(2185)
# # [1] 0.50.8
# message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")
# 
# p2<-ggplot(data=dat2_outlier, aes(V5)) + 
#   geom_histogram(bins = 40, color="black", fill="white")+
#   scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~wild~transect~"("~ARN~-~COH~")")) + 
#   #ggtitle("Same SNPs in ARN - COH but ratio \nis 1109 out of 2185 (0.51) ")+
#   ggtitle("Same SNPs in ARN - COH but ratio \nis 299 out of 540 (0.55) ")+
#   theme(axis.text=element_text(size=14),
#         text = element_text(size=14), # add family="times" in macOS for time font
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.5),
#         plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname3 = "ps_Del20_challenge.txt"
dat3 = read.delim(pname3, header = FALSE, sep='\t')
dat3$id1 = paste0(dat3$V1,'_',dat3$V2)
dat3_outlier = dat3[which(dat3$id1 %in% outlier_id),]
neg <- length(dat3_outlier$V5[which(dat3_outlier$V5<0)])
pos <- length(dat3_outlier$V5[which(dat3_outlier$V5>0)])
# 1092/(2185)
# [1] 0.50
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p3<-ggplot(data=dat3_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~"2020"~challenge)) + 
  ggtitle("Same SNPs in 2020 challenge but ratio \nis 1092 out of 2185 (0.50) ")+
  #ggtitle("Same SNPs in 2020 challenge but ratio \nis 273 out of 540 (0.51) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

ggarrange(p, p1, p2, p3,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

dev.off()

###########################################################################
## Check SNP direction in independent SGS using wild transect as a prior ## 
########################################################################### 

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/indep_SGS_ps")

# check the delta_p directionality
#check_direct <- function(pname){
pname = "ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
length(dat$V5[which(dat$V5<0)])
length(dat$V5[which(dat$V5>=0)])
adj = p.adjust(dat$V6, method = 'BH')
length(dat$V5[which(adj< 0.05)])
dat_outlier <- dat[which(adj< 0.05),]
outlier_id <- paste0(dat_outlier$V1,'_',dat_outlier$V2)
neg <- length(dat_outlier$V5[which(dat_outlier$V5<0)])
pos <- length(dat_outlier$V5[which(dat_outlier$V5>=0)])
# 1805/(1805+380)
# [1] 0.826087
#share_all = sum(sign(as.numeric(dp1)) == sign(as.numeric(dp2)))
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")
#}
par(mar=c(4,4,2,8))
par(mfrow=c(2,1))
jpeg("SNP_FDR05_Del19.jpg", width = 12, height = 8, units = 'in', res = 300)

p<-ggplot(data=dat_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(SNP~candidates~Delta~italic(p)~"in"~wild~transect~"("~HC~-~NB~")")) + 
  ggtitle("Candidates with positive changes is \n1767 out of 2073 (0.85)")+
  #ggtitle("Candidates with positive changes is \n441 out of 540 (0.82)")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname1 = "ps_Del19_challenge.txt"
dat1 = read.delim(pname1, header = FALSE, sep='\t')
dat1$id1 = paste0(dat1$V1,'_',dat1$V2)
dat1_outlier = dat1[which(dat1$id1 %in% outlier_id),]
neg <- length(dat1_outlier$V5[which(dat1_outlier$V5<0)])
pos <- length(dat1_outlier$V5[which(dat1_outlier$V5>0)])
# 1042/(2185)
# [1] 0.4768
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p1<-ggplot(data=dat1_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~'2019'~challenge)) + 
  ggtitle("Same SNPs in 2019 challenge but ratio \nis 1041 out of 2073 (0.50) ")+
  #ggtitle("Same SNPs in HC - NB but ratio \nis 248 out of 540 (0.46) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname2 = "ps_Del19_SR_NB.txt"
dat2 = read.delim(pname2, header = FALSE, sep='\t')
dat2$id1 = paste0(dat2$V1,'_',dat2$V2)
dat2_outlier = dat2[which(dat2$id1 %in% outlier_id),]
neg <- length(dat2_outlier$V5[which(dat2_outlier$V5<0)])
pos <- length(dat2_outlier$V5[which(dat2_outlier$V5>0)])
# 1109/(2185)
# [1] 0.50.8
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p2<-ggplot(data=dat2_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~wild~transect~"("~SR~-~NB~")")) + 
  ggtitle("Same SNPs in SR - NB but ratio \nis 1733 out of 2073 (0.84) ")+
  #ggtitle("Same SNPs in SR - NB but ratio \nis 299 out of 540 (0.55) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

# pname2 = "ps_Del19_ARN_COH.txt"
# dat2 = read.delim(pname2, header = FALSE, sep='\t')
# dat2$id1 = paste0(dat2$V1,'_',dat2$V2)
# dat2_outlier = dat2[which(dat2$id1 %in% outlier_id),]
# neg <- length(dat2_outlier$V5[which(dat2_outlier$V5<0)])
# pos <- length(dat2_outlier$V5[which(dat2_outlier$V5>0)])
# # 1109/(2185)
# # [1] 0.50.8
# message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")
# 
# p2<-ggplot(data=dat2_outlier, aes(V5)) + 
#   geom_histogram(bins = 40, color="black", fill="white")+
#   scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~wild~transect~"("~ARN~-~COH~")")) + 
#   #ggtitle("Same SNPs in ARN - COH but ratio \nis 1109 out of 2185 (0.51) ")+
#   ggtitle("Same SNPs in ARN - COH but ratio \nis 299 out of 540 (0.55) ")+
#   theme(axis.text=element_text(size=14),
#         text = element_text(size=14), # add family="times" in macOS for time font
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.5),
#         plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

pname3 = "ps_Del20_challenge.txt"
dat3 = read.delim(pname3, header = FALSE, sep='\t')
dat3$id1 = paste0(dat3$V1,'_',dat3$V2)
dat3_outlier = dat3[which(dat3$id1 %in% outlier_id),]
neg <- length(dat3_outlier$V5[which(dat3_outlier$V5<0)])
pos <- length(dat3_outlier$V5[which(dat3_outlier$V5>0)])
# 1092/(2185)
# [1] 0.50
message("number of SNPs with positive delta_p is ", pos, " out of ", (pos+neg), " (", pos/(pos+neg), ")")

p3<-ggplot(data=dat3_outlier, aes(V5)) + 
  geom_histogram(bins = 40, color="black", fill="white")+
  scale_x_continuous(expression(Same~SNPs~Delta~italic(p)~"in"~"2020"~challenge)) + 
  ggtitle("Same SNPs in 2020 challenge but ratio \nis 1033 out of 2073 (0.50) ")+
  #ggtitle("Same SNPs in 2020 challenge but ratio \nis 273 out of 540 (0.51) ")+
  theme(axis.text=element_text(size=14),
        text = element_text(size=14), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm")) 

ggarrange(p, p1, p2, p3,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

dev.off()

#########################################
## Outlier SNPs test for directionality #
######################################### 

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/p_combine/SR_REF19_ARN_COH")
# check the delta_p directionality

check_direct <- function(pname){
  #pname = "./REF19_CHR19_NB_HC_out_all_fish.txt"
  dat = read.delim(pname, header = TRUE, sep='\t')
  dp1 = dat$dp_c
  dp2 = dat$dp_w
  share_all = sum(sign(as.numeric(dp1)) == sign(as.numeric(dp2)))
  message("number of SNPs with same directionality is ", share_all, " out of ", dim(dat)[1], " (", share_all/dim(dat)[1], ")")
}
check_direct("SR_REF19_ARN_COH_out_all_fish.txt")
check_direct("SR_REF19_ARN_COH_out_0.05_fish.txt")
check_direct("REF19_CHR19_ARN_COH_out_all_fish.txt")
check_direct("REF19_CHR19_ARN_COH_out_0.05_fish.txt")
check_direct("REF19_CHR19_NB_HC_out_all_fish.txt")
check_direct("REF19_CHR19_NB_HC_out_0.05_fish.txt")
check_direct("REF19_CHR19_NB_HC_out_0.01_fish.txt")
check_direct("REF19_CHR19_REF20_CHR20_out_0.05_fish.txt")
check_direct("HC_CS19_CLP_HC-VA_out_0.05_fish.txt")

cnt_total = as.data.frame(seq(1,2335739,1))
outlier = 6803

library(dplyr)
n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
pname = "./REF19_CHR19_NB_HC_out_all_fish.txt"
dat = read.delim(pname, header = TRUE, sep='\t')
for (i in 1:n_bootstraps){
  
  idx = sample_n(cnt_total, 6803)
  dat_ = dat[idx$`seq(1, 2335739, 1)`,]
  dp1_ = dat_$dp_c
  dp2_ = dat_$dp_w
  cnt = sum(sign(as.numeric(dp1_)) == sign(as.numeric(dp2_)))
  boot_cnt[i] = cnt
}
boot_cnt
p = sum(boot_cnt >= 4608)/length(boot_cnt)
p
hist(boot_cnt)
abline(v = mean(boot_cnt), col = "blue", lwd = 2)
abline(v = 4608, col = "red", lwd = 2, lty=2)
# This test confirms whether the normal distribution of the data is violated.
shapiro.test(boot_cnt)
t.test(boot_cnt,conf.level=0.95)

#####################################
##########  2d density plot  ########
#####################################
library(ggplot2)
library(showtext)
library(MASS)
library(tidyverse)
library(ggpubr)

font_add_google("Noto Sans SC")
showtext_auto()

# highlight1 = read.delim("GEA_BF_20_Del19_FDR_2K.intersect", header = FALSE, sep='\t')
# highlight_plot1 = outlier_list[which(outlier_list$id %in% highlight1$V4),] # note the snp count is less than GEA_BF_20_Del19_FDR_10K.intersect because of duplicate intersects
# 
# SGS_Del19_FDR<-read.table("Del19_FDR_outlier.list", header = T)
# GEA_BF20<-read.table("salinity_2032113_outlier_BF20.txt", header = T)
# GEA_BF20$id = paste0(GEA_BF20$chromo,'_',GEA_BF20$position)
# SGS_Del19_FDR$id[which(SGS_Del19_FDR$id %in% GEA_BF20$id)]
# 
# highlight2 = read.delim("GEA_BF_20_Del19_FDR_2K.intersect", header = FALSE, sep='\t')
# highlight_plot2 = outlier_list[which(outlier_list$id %in% "NC_035786.1_8772371"),] # note the snp count is less than GEA_BF_20_Del19_FDR_10K.intersect because of duplicate intersects


pname = "./REF19_CHR19_ARN_COH_out_0.05_fish.txt"
dat = read.delim(pname, header = TRUE, sep='\t')

highlight = "./REF19_CHR19_NB_HC_out_0.01_fish.txt"
highlight_plot1 = read.delim(highlight, header = TRUE, sep='\t')
# https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2.html#hex
# Area + contour
tiff("11_SGS_indep_p_zcomb_FDR005.tiff", units="in", width=8, height=8, res=300)
p1 <- ggplot(dat, aes(x=dp_c, y=dp_w) ) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme(legend.position='right')+
  geom_point(color="dimgray", alpha = 0.2)+
  geom_density_2d(colour="white", size=0.8, alpha = 0.5) +
  #stat_density2d(aes(color = ..level..))+
  #stat_density_2d(aes(fill = factor(stat(level))), geom = "polygon",alpha = 0.5, colour="grey") +
  #geom_point(data=highlight_plot1, aes(x=dp_c, y=dp_w) , color='red', size=2) +
  #geom_point(data=highlight_plot2, aes(x=Cdelta_p, y=Wdelta_p) , color='blue', size=2) +
  scale_x_continuous(expression(Delta~italic(p)~"in"~'2019'~challenge~"("~Surv~-~Ref~")"), limits = c(-0.5, 0.5)) + 
  scale_y_continuous(expression(Delta~italic(p)~"in"~wild~transect~"("~ARN~-~COH~")"), limits = c(-0.5, 0.5))+ 
  
  #scale_x_continuous(expression(Delta~italic(p)~"in"~challenge~"("~Surv~-~Ref~")"), limits = c(-0.5, 0.5)) + 
  #scale_y_continuous(expression(Delta~italic(p)~"in"~wild~transect~"("~HC~-~NB~")"), limits = c(-0.5, 0.5))+ 
  
  #scale_x_continuous(expression(Delta~italic(p)~"in"~'2019'~challenge~"("~Surv~-~Ref~")"), limits = c(-0.5, 0.5)) + 
  #scale_y_continuous(expression(Delta~italic(p)~"in"~'2020'~challenge~"("~Surv~-~Ref~")"), limits = c(-0.5, 0.5))+
  
  theme_pubr() + 
  theme(axis.text=element_text(size=20),
        text = element_text(size=20), # add family="times" in macOS for time font
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))
p1 
#p1 + stat_density_2d_filled(alpha = 0.3, contour_var="density") 
dev.off()

###############################################
##########  delta_p distribution plot  ########
###############################################
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(grid)
library(gtable)
library(export)
pname = "./REF19_CHR19_NB_HC_out_0.01_fish.txt"
pname = "./REF19_CHR19_ARN_COH_out_0.05_fish.txt"
dat = read.delim(pname, header = TRUE, sep='\t')
dp1 = dat$dp_c
dp2 = dat$dp_w
share_all = sign(as.numeric(dp1)) == sign(as.numeric(dp2))
dat_same = dat[share_all,]
dat_same$group = "outliers with same directionality"
write.table(dat_same, "./REF19_CHR19_ARN_COH_FDR0.05_same.txt", row.names=F, quote=F, sep="\t")

pname_all = "REF19_CHR19_ARN_COH_out_all_fish.txt"
dat_all = read.delim(pname_all, header = TRUE, sep='\t')
dat_all$group = "genome-wide"

dat_plot = rbind(dat_same, dat_all)

p2 <- ggplot(dat_plot, aes(x=dp_c, fill=group)) +  
  geom_density(alpha = 0.4) + theme_void() + theme(legend.position = "right") + scale_fill_manual(values=c("#F5F5F5", "#666666"))
  
p2 

p3 <- ggplot(dat_plot, aes(x=dp_w, fill=group)) +  
  geom_density(alpha = 0.4) + theme_void() + theme(legend.position = "none") + scale_fill_manual(values=c("#F5F5F5", "#666666")) + coord_flip()

p3

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
tiff("11_SGS_indep_p_zcomb_FDR005_1.tiff", units="in", width=8, height=8, res=300)
comb_plot<-grid.arrange(p2,empty,p1, p3, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()
gtable_show_layout(comb_plot)

graph2ppt(file="challenge_plot_FDR01.pptx", width=10, height=10)
annotate_figure(figure,
                top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
                bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                right = "I'm done, thanks :-)!",
                fig.lab = "Figure 1", fig.lab.face = "bold"
)



p <- qplot(1,1)
g <- ggplotGrob(p)

panel_id <- g$layout[g$layout$name == "panel",c("t","l")]
g <- gtable_add_cols(g, unit(1,"cm"))

g <- gtable_add_grob(g, rectGrob(gp=gpar(fill="red")),
                     t = panel_id$t, l = ncol(g))

g <- gtable_add_rows(g, unit(1,"in"), 0)
g <- gtable_add_grob(g, rectGrob(gp=gpar(fill="blue")),
                     t = 1, l = panel_id$l)

grid.newpage()
grid.draw(g)


##################################################
##########  starting p distribution plot  ########
###################################
###############
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(grid)
library(gtable)
library(export)
pname = "./ps_NB.txt"
dat = read.delim(pname, header = TRUE, sep='\t')
head(dat)
p2 <- ggplot(dat, aes(x=p_NB, color=group)) +  
  geom_density(alpha = 0.2, size = 1.2) + theme(legend.position = "right") +
  #geom_vline(aes(xintercept=grp.mean, color=group), linetype="dashed")
  scale_color_manual(values=c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600")) +
  ylab("Density") + xlab("Starting allele frequency in 2019 New Bed (NB)") +
  theme_classic()

p2 

