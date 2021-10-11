####################################
##########  log(FDR) plot ##########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/20_Balancing_selection/plot/chr5/mafs")
source("manhattan.R")

jpeg("Mahattan_bs_SGS.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mar=c(4,8,1,6))
par(mfrow=c(3,1))

# B0 scores
chr = "10"
file = 'NB_NC_035789.1_dafs_1kb_500b.txt' 
bs = read.delim(file, header = TRUE, sep='\t')
bs$chromo = chr
colnames(bs)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat', 'nSite', 'chromo')
bs$SNP = paste0(bs$chromo,'_',bs$position)
bs$chromo = as.numeric(bs$chromo)

# B0MAF scores
chr = "10"
file = 'NB_NC_035789.1.mafs.output.1kb_500b.txt' 
bs_maf = read.delim(file, header = TRUE, sep='\t')
bs_maf$chromo = chr
colnames(bs_maf)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat', 'nSite', 'chromo')
bs_maf$SNP = paste0(bs_maf$chromo,'_',bs_maf$position)
bs_maf$chromo = as.numeric(bs_maf$chromo)

name1 = 'outlier_pos.txt' 
h1 = read.delim(name1, header = FALSE, sep=',')
h1 = as.list(h1)

name2 = 'outlier_neg.txt' 
h2 = read.delim(name2, header = FALSE, sep=',')
h2 = as.list(h2)

file = 'SGS_REF19_NB_chr10.txt' 
dat = read.delim(file, header = TRUE, sep='\t')
dat$logFDR = -log(dat$FDR) 
dat$SNP <- paste0(dat$chromo,'_',dat$position)

manhattan(chr="chromo",bp="position",p="CLR", snp = "SNP", bs, logp=FALSE, cex.axis = 1.2, ylim = c(0, max(bs$CLR)),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="B0 score", cex.lab=1.5) 

manhattan(chr="chromo",bp="position",p="CLR", snp = "SNP", bs_maf, logp=FALSE, cex.axis = 1.2, ylim = c(0, max(bs$CLR)),
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="B0maf score", cex.lab=1.5) 

manhattan(chr="chromo",bp="position",p="logFDR", snp = "SNP", dat, highlight1 = h1$V1, highlight2 = h2$V1, logp=FALSE, cex.axis = 1.2,
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab=expression(log~FDR), cex.lab=1.5) 

dev.off()


#######################
#  plot delta_p      #
#######################
# plot distribution

#jpeg("EA_60_outlier_27_direction_trend.jpg", width = 16, height = 9, units = 'in', res = 300)
#par(mfrow=c(1,1))
outlier <- read.delim("/Users/ryan/Documents/Ryan_workplace/DelBay_adult/08_fish_exact/plot_allele_frequency_trend/REF19_CHR19_NB_HC_out_0.05_fish_new.txt", header = FALSE, sep='\t')
outlier$id = paste0(outlier$V1,'_',outlier$V2)

# load the dataset
# load the snp info and allele frequency data from each population 
fname = "by_pop_0.05_pctind0.7_maxdepth3.mafs"  
df_maf = read.delim(fname, header = TRUE, sep=' ')
head(df_maf)
dim(df_maf)
df_maf$id = paste0(df_maf$chromo,'_',df_maf$position)
df_maf <- df_maf[with(df_maf, order(id)),]
head(df_maf)
# load the delta_p information from DelBay19 challenge
fname = "ps_Del19_challenge.txt"  
df_deltap = read.delim(fname, header = FALSE, sep='\t')
dim(df_deltap)
colnames(df_deltap) <- c("chromo", "position", "p1", "p0", "delta_p", "ps", "outlier")
df_deltap$id = paste0(df_deltap$chromo,'_',df_deltap$position)
df_deltap <- df_deltap[with(df_deltap, order(id)),]
head(df_deltap)
sum(df_maf$id == df_deltap$id)

data1 = cbind(df_maf$chromo[df_maf$id %in% outlier$id], 
              df_maf$position[df_maf$id %in% outlier$id], 
              df_maf$major[df_maf$id %in% outlier$id], 
              df_maf$minor[df_maf$id %in% outlier$id],
              df_maf$id[df_maf$id %in% outlier$id],
              df_maf$freqHC[df_maf$id %in% outlier$id],
              df_maf$freqARN[df_maf$id %in% outlier$id],
              df_maf$freqCOH[df_maf$id %in% outlier$id],
              df_maf$freqSR[df_maf$id %in% outlier$id],
              df_maf$freqNB[df_maf$id %in% outlier$id],
              df_maf$freqHC[df_maf$id %in% outlier$id]-df_maf$freqNB[df_maf$id %in% outlier$id])

data2 =  cbind(df_deltap$id[df_deltap$id %in% outlier$id],
               df_deltap$p1[df_deltap$id %in% outlier$id],
               df_deltap$p0[df_deltap$id %in% outlier$id],
               df_deltap$delta_p[df_deltap$id %in% outlier$id],
               df_deltap$ps[df_deltap$id %in% outlier$id])

outlier_list <- as.data.frame(cbind(data1, data2))
# rename the columns
colnames(outlier_list) <- c("chromo", "position", "major", "minor", "id", 
                            "freqHC", "freqARN", "freqCOH", "freqSR", "freqNB", "Wdelta_p",
                            "id2", "freqSurv", "freqRef", "Cdelta_p", "ps" )
# change the column format
i <- c(seq(6,11), seq(13,16))
outlier_list[ , i] <- apply(outlier_list[ , i], 2,            # Specify own function within apply
                            function(x) as.numeric(as.character(x)))
sapply(outlier_list, class)
# take a look at outlier_list
head(outlier_list)
# polarize the "minor" alleles
for (i in 1:length(outlier_list$position)){
  if (outlier_list$Cdelta_p[i] < 0){
    tmp = outlier_list$major[i]
    outlier_list$major[i] = outlier_list$minor[i]
    outlier_list$minor[i] = tmp
    outlier_list$freqHC[i] = 1-outlier_list$freqHC[i]
    outlier_list$freqARN[i] = 1-outlier_list$freqARN[i]
    outlier_list$freqCOH[i] = 1-outlier_list$freqCOH[i]
    outlier_list$freqSR[i] = 1-outlier_list$freqSR[i]
    outlier_list$freqNB[i] = 1-outlier_list$freqNB[i]
    outlier_list$Wdelta_p[i] = -outlier_list$Wdelta_p[i]
    outlier_list$freqSurv[i] = 1-outlier_list$freqSurv[i]
    outlier_list$freqRef[i] = 1-outlier_list$freqRef[i]
    outlier_list$Cdelta_p[i] = -outlier_list$Cdelta_p[i]
  }
}
# take a look at outlier_list after polarization
head(outlier_list)
df_same_direction = outlier_list
# rank order based NB frequency
df_same_direction_rank <- df_same_direction[with(df_same_direction, order(freqRef)),]
# customize the color
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}
mycol <- t_col("black", perc = 50, name = "lt.grey")

hist(df_same_direction_rank$Wdelta_p, main=NULL, breaks = 50,
     xlab = expression(Delta~italic(p)), ylab = "Frequency")
# plot the trend
cnt = length(df_same_direction_rank$chromo)
cnt

jpeg("Fisher_18_outlier_direction_trend.tiff", width = 16, height = 9, units = 'in', res = 300)
#par(mfrow=c(1,1))

order = seq(1,length(df_same_direction_rank$chromo),1)
df_same_direction_rank$order = order
plot(df_same_direction_rank$order, df_same_direction_rank$freqNB,col="red",pch=20,cex=2.2,ylim=c(0,1),xlab="Outliers",ylab="Allele frequency",xlim=c(0,cnt)) # ylab=expression(italic("p"))
for (l in 1:nrow(df_same_direction_rank))
{
  segments(l,df_same_direction_rank$freqHC[l],l,df_same_direction_rank$freqNB[l],col=mycol,lwd=2.0, lty = "dotted")
  if (df_same_direction_rank$Wdelta_p[l]<=0)
  {
    points(l,df_same_direction_rank$freqHC[l],col="blue",pch=6,cex=2.2)
  }
  if (df_same_direction_rank$Wdelta_p[l]>0)
  {
    points(l,df_same_direction_rank$freqHC[l],col="blue",pch=17,cex=2.2)
  }
}

points(df_same_direction_rank$order+0.24, df_same_direction_rank$freqRef,col="red",pch=20,cex=2.2,ylim=c(0,1),ylab=expression(italic("p")),xlab="SNPs",xlim=c(0,cnt))
for (l in 1:nrow(df_same_direction_rank))
{
  segments(l+0.24,df_same_direction_rank$freqSurv[l],l+0.24,df_same_direction_rank$freqRef[l],col=mycol,lwd=2.0, lty = "dotted") #@ add by HG
  if (df_same_direction_rank$Wdelta_p[l]<=0)
  {
    points(l+0.24,df_same_direction_rank$freqSurv[l],col="green",pch=6,cex=2.2)
  }
  if (df_same_direction_rank$Wdelta_p[l]>0)
  {
    points(l+0.24,df_same_direction_rank$freqSurv[l],col="green",pch=17,cex=2.2) #@ add by HG
  }
}

dev.off()


