
#######################
#  Adjust p-value    #
#######################

####################################
#  extract p<0.05 and FDR < 0.1    #
####################################

ps_adj <- function(pname, outname1, outname2) {
  #pname = "ps_18_SGS_HC_NB.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  # replace chromosome if it is numerical
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  # count how many SNPs with p < 0.05
  ps = dat$V6
  ps_cut <- 1e-4
  ps_cnt = length(ps[ps<ps_cut])
  message(paste0("p-value < ", ps_cut, ": "), ps_cnt)
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$V1,'_',dat$V2)
  ps_outlier <- dat[which(dat$V6<ps_cut),][,1:5]
  pos_cnt = length(ps_outlier[which(ps_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(ps_outlier[which(ps_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  # adjust the p-values with FDR correction
  dat$adj = p.adjust(ps, method = 'BH')
  adj=dat$adj
  # output a table for SNPs with p < 0.05
  outlier1 <- as.data.frame(cbind(dat$V1[ps<ps_cut], dat$V2[ps<ps_cut], dat$V3[ps<ps_cut], dat$V4[ps<ps_cut], dat$V5[ps<ps_cut], dat$V6[ps<ps_cut], dat$adj[ps<ps_cut], dat$id[ps<ps_cut]))
  colnames(outlier1) <- c("chr","pos","p1","p2","deltap","p_value","fdr", "id")
  outlier1 <- outlier1[with(outlier1, order(chr, pos)),]
  write.table(outlier1,outname1, row.names = FALSE, sep="\t", quote = FALSE)
  # output a table for SNPs with FDR < 0.05
  adj_cut <- 0.1
  outlier2 <- as.data.frame(cbind(dat$V1[adj<adj_cut], dat$V2[adj<adj_cut], dat$V3[adj<adj_cut], dat$V4[adj<adj_cut], dat$V5[adj<adj_cut], dat$V6[adj<adj_cut], dat$adj[adj<adj_cut], dat$id[adj<adj_cut]))
  #outlier2 <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], dat$V4[adj<0.05], dat$V5[adj<0.05], dat$V6[adj<0.05], dat$adj[adj<0.05], dat$id[adj<0.05]))
  colnames(outlier2) <- c("chr","pos","p1","p2","deltap","p_value","fdr", "id")
  outlier2 <- outlier2[with(outlier2, order(chr, pos)),]
  write.table(outlier2,outname2, row.names = FALSE, sep="\t", quote = FALSE)
  # count how many SNPs with FDR < adj_cut
  message(paste0("FDR < ", adj_cut, ": ", length(outlier2$id)))
  # count how many SNPs with positive vs negative frequency change
  adj_outlier <- dat[which(dat$adj<adj_cut),][,1:5]
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
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/11_SGS/two_sided_p_value_testassoc_list/")
ps_adj("ps_18_SGS_HC_NB.txt", "18_SGS_HC_NB_ps_outlier.list", "18_SGS_HC_NB_FDR_outlier.list")
ps_adj("ps_19_SGS_HC_NB.txt", "19_SGS_HC_NB_ps_outlier.list", "19_SGS_HC_NB_FDR_outlier.list")
ps_adj("ps_21_SGS_HC_NB.txt", "21_SGS_HC_NB_ps_outlier.list", "21_SGS_HC_NB_FDR_outlier.list")
ps_adj("ps_19_SGS_Sur_Ref.txt", "19_SGS_Sur_Ref_ps_outlier.list", "19_SGS_Sur_Ref_FDR_outlier.list")
ps_adj("ps_20_SGS_Sur_Ref.txt", "20_SGS_Sur_Ref_ps_outlier.list", "20_SGS_Sur_Ref_FDR_outlier.list")
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

SGS_HC_NB_18 = make_ID('18_SGS_HC_NB_FDR_outlier.list')
SGS_HC_NB_19 = make_ID('19_SGS_HC_NB_FDR_outlier.list')
SGS_HC_NB_19_flip = make_ID('19_SGS_HC_NB_FDR_outlier_flip.list')
SGS_HC_NB_21 = make_ID('21_SGS_HC_NB_FDR_outlier.list')
SGS_Sur_Ref_19 = make_ID('19_SGS_Sur_Ref_FDR_outlier.list')
SGS_Sur_Ref_20 = make_ID('20_SGS_Sur_Ref_FDR_outlier.list')


SGS_HC_NB_18 = make_ID('18_SGS_HC_NB_ps_outlier.list')
SGS_HC_NB_19 = make_ID('19_SGS_HC_NB_ps_outlier.list')
#SGS_HC_NB_19_flip = make_ID('19_SGS_HC_NB_ps_outlier_flip.list')
SGS_HC_NB_21 = make_ID('21_SGS_HC_NB_ps_outlier.list')
SGS_Sur_Ref_19 = make_ID('19_SGS_Sur_Ref_ps_outlier.list')
SGS_Sur_Ref_20 = make_ID('20_SGS_Sur_Ref_ps_outlier.list')


x <- list(
  A = SGS_HC_NB_18, 
  B = SGS_HC_NB_19, 
  C = SGS_HC_NB_21,
  D = SGS_Sur_Ref_19,
  E = SGS_Sur_Ref_20
)

library("ggvenn")
names(x) <- c("HC_NB_18","HC_NB_19","HC_NB_21", "Sur_Ref_19", "Sur_Ref_20")
fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#bc5090")

ggvenn(
  x, columns = c("HC_NB_18", "HC_NB_19", "HC_NB_21", "Sur_Ref_19"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggvenn(
  x, columns = c("HC_NB_19", "Sur_Ref_19"),
  fill_color = c("#EFC000FF", "#bc5090"),
  stroke_size = 0.5, set_name_size = 4
)

ggvenn(
  x, columns = c("HC_NB_18", "HC_NB_19", "HC_NB_21"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5
)

ggvenn(
  x, columns = c("Sur_Ref_19", "Sur_Ref_20"),
  fill_color = c("#CD534CFF", "#bc5090"),
  stroke_size = 0.5
)

challenge_shared<- intersect(SGS_Sur_Ref_19, SGS_Sur_Ref_20)
write.table(challenge_shared,"challenge_shared.bed", row.names = FALSE, sep="\t", quote = FALSE)

intersect(SGS_HC_NB_21, intersect(SGS_HC_NB_18, SGS_HC_NB_19))

length(SGS_HC_NB_19)
length(SGS_HC_NB_19_flip)
length(intersect(SGS_HC_NB_19_flip, SGS_HC_NB_19))
# not so many overlaps given because it has different assumption
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





