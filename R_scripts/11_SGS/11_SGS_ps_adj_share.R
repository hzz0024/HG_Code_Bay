setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/11_SGS")
#######################
#  Adjust p-value    #
#######################

# load the p-values (2032113 SNPs, from global shared 1x coverage data)

####################################
#  extract p<0.05 and FDR < 0.05   #
####################################

ps_adj <- function(pname, outname1, outname2) {
  #pname = "ps_Del19_challenge.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  # replace chromosome if it is numerical
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  # count how many SNPs with p < 0.05
  ps = dat$V6
  ps_cnt = length(ps[ps<0.05])
  message("p-value < 0.05:", ps_cnt)
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$V1,'_',dat$V2)
  ps_outlier <- dat[which(dat$V6<0.05),][,1:5]
  pos_cnt = length(ps_outlier[which(ps_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(ps_outlier[which(ps_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  # adjust the p-values with FDR correction
  dat$adj = p.adjust(ps, method = 'BH')
  adj=dat$adj
  # output a table for SNPs with p < 0.05
  outlier1 <- as.data.frame(cbind(dat$V1[ps<0.05], dat$V2[ps<0.05], dat$V3[ps<0.05], dat$V4[ps<0.05], dat$V5[ps<0.05], dat$V6[ps<0.05], dat$adj[ps<0.05], dat$id[ps<0.05]))
  colnames(outlier1) <- c("chr","pos","p1","p2","deltap","p_value","fdr", "id")
  outlier1 <- outlier1[with(outlier1, order(chr, pos)),]
  write.table(outlier1,outname1, row.names = FALSE, sep="\t", quote = FALSE)
  # output a table for SNPs with FDR < 0.05
  outlier2 <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], dat$V4[adj<0.05], dat$V5[adj<0.05], dat$V6[adj<0.05], dat$adj[adj<0.05], dat$id[adj<0.05]))
  colnames(outlier2) <- c("chr","pos","p1","p2","deltap","p_value","fdr", "id")
  outlier2 <- outlier2[with(outlier2, order(chr, pos)),]
  write.table(outlier2,outname2, row.names = FALSE, sep="\t", quote = FALSE)
  # count how many SNPs with FDR < 0.05
  message("FDR < 0.05:" ,length(outlier2$id))
  # count how many SNPs with positive vs negative frequency change
  adj_outlier <- dat[which(dat$adj<0.05),][,1:5]
  pos_cnt = length(adj_outlier[which(adj_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(adj_outlier[which(adj_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  # depict the p-value and FDR distribution
  jpeg(paste0(pname, ".jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(1,2))
  hist(ps, main="", xlab = "p-value", ylab = "Frequency")
  hist(adj, main="", xlab = "FDR", ylab = "Frequency")
  dev.off()
}

ps_adj("ps_Del19_challenge.txt", "challenge_ps_outlier.list", "challenge_FDR_outlier.list")
ps_adj("ps_Del19_HC_NB.txt", "HC_NB_ps_outlier.list", "HC_NB_FDR_outlier.list")
ps_adj("ps_Del19_HC_SR.txt", "HC_SR_ps_outlier.list", "HC_SR_FDR_outlier.list")


#######################
#  Shared outliers    #
#######################

make_ID <- function(file_name){
  dat = read.delim(file_name, header = TRUE, sep='\t')
  dat = dat[with(dat, order(chr, pos)),]
  dat_ID = dat$id
  message(length(dat$id))
  return(dat_ID)
}

Del19_p5 = make_ID('challenge_ps_outlier.list')
HC_NB_p5 = make_ID('HC_NB_ps_outlier.list')
HC_SR_p5 = make_ID('HC_SR_ps_outlier.list')

a <- intersect(Del19_p5, HC_SR_p5)
length(a)
b <- intersect(Del19_p5, HC_NB_p5)
length(b)

i <- intersect(HC_SR_p5, (intersect(Del19_p5, HC_NB_p5)))
length(i)

Del19 = make_ID('challenge_FDR_outlier.list')
HC_NB = make_ID('HC_NB_FDR_outlier.list')
HC_SR = make_ID('HC_SR_FDR_outlier.list')

a <- intersect(Del19, HC_SR)
length(a)
b <- intersect(Del19, HC_NB)
length(b)

i <- intersect(HC_SR, (intersect(Del19, HC_NB)))
length(i)



#######################
#  Check delta_p      #
#######################

# check outlier delta_p in each contrast

check_outlier <- function(pname){
  #pname = "ps_Del19_challenge.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$V1,'_',dat$V2)
  outlier <- dat[which(dat$id %in% c),][,1:5] # modify the i to other common shared outliers
  colnames(outlier) <- c("chr","pos","p1","p2","deltap")
  cnt = length(outlier$deltap)
  pos_cnt = length(outlier[which(outlier$deltap>0),]$pos)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(outlier[which(outlier$deltap<0),]$pos)
  message("negative change SNPs: ", neg_cnt)
  #outlier1 <- outlier[with(outlier, order(-deltap)),] # order by deltap from high to low
  outlier1 <- outlier[with(outlier, order(chr, pos)),] # order by chromosome and position
  return(outlier1)
}

check_outlier("ps_Del19_challenge.txt")
check_outlier("ps_Del19_HC_NB.txt")
check_outlier("ps_Del19_HC_SR.txt")

pname = "HC_NB_FDR_outlier.list" # 2-side p-value
dat = read.delim(pname, header = TRUE, sep='\t')
cnt = length(dat$chr)
cnt
outlier=dat[,1:5]
colnames(outlier) <- c("chr","pos","p1","p2","deltap")
pos_cnt = length(outlier[which(outlier$deltap>0),]$pos)
message("positive change SNPs: ", pos_cnt)
neg_cnt = length(outlier[which(outlier$deltap<0),]$pos)
message("negative change SNPs: ", neg_cnt)
outlier <- outlier[with(outlier, order(-deltap)),]
#######################
#  plot delta_p      #
#######################
# plot distribution
jpeg("SGS_HC_SR_ps_79542.jpg", width = 16, height = 9, units = 'in', res = 300)
jpeg("SGS_shared_823.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,2))

hist(outlier$deltap, main=NULL, breaks = 50,
     xlab = expression(Delta~italic(p)), ylab = "Frequency")
# plot the trend
plot(outlier$p2,col="red",pch=".",cex=3,ylim=c(0,1),ylab=expression(italic("p")),xlab="SNPs",xlim=c(0,cnt+1))
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
mycol <- t_col("grey", perc = 50, name = "lt.grey")
for (l in 1:nrow(outlier))
{
  segments(l,outlier$p1[l],l,outlier$p2[l],col=mycol,lwd=1.2, lty = "dotted")
  if (outlier$deltap[l]<=0)
  {
    points(l,outlier$p1[l],col="green",pch=6,cex=0.7)
  }
  if (outlier$deltap[l]>0)
  {
    points(l,outlier$p1[l],col="blue",pch=2,cex=0.7)
  }
}
abline(v=pos_cnt+0.5,lty=2, col="darkgrey")
dev.off()





