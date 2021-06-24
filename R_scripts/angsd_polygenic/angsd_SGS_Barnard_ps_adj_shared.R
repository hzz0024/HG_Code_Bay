#######################
##  Adjust p-value   ##
#######################
setwd("~/Documents/Ryan_workplace/DelBay_adult/19_share_SGS_barnard/")
ps_adj <- function(pname, outname1, outname2) {
  #pname = "barnard_REF19_CHR19.txt" # 2-side p-value
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

ps_adj("REF19_CHR19_out_all_z.txt", "Del19_ps_outlier.list", "Del19_FDR_outlier.list")
ps_adj("REF20_CHR20_out_all_z.txt", "Del20_ps_outlier.list", "Del20_FDR_outlier.list")
ps_adj("SR_HC_out_all_z.txt", "SR_HC_ps_outlier.list", "SR_HC_FDR_outlier.list")
ps_adj("NB_HC_out_all_z.txt", "NB_HC_ps_outlier.list", "NB_HC_FDR_outlier.list")

#######################
#  Shared outliers    #
#######################

make_ID <- function(file_name){
  #file_name ="Del19_outlier.list"
  dat = read.delim(file_name, header = TRUE, sep='\t')
  dat = dat[with(dat, order(chr, pos)),]
  message(length(dat$id))
  return(dat$id)
}

Del19 = make_ID('Del19_FDR_outlier.list')
Del20 = make_ID('Del20_FDR_outlier.list')
HC_SR = make_ID('SR_HC_FDR_outlier.list')
HC_NB = make_ID('NB_HC_FDR_outlier.list')

a <- intersect(Del19, Del20)
length(a)
b <- intersect(Del19, HC_SR)
length(b)
c <- intersect(Del19, HC_NB)
length(c)
e <- intersect(Del20, HC_SR)
length(e)
d <- intersect(Del20, HC_NB)
length(d)
f <- intersect(HC_NB, HC_SR)
length(f)

Del19 = make_ID('Del19_ps_outlier.list')
Del20 = make_ID('Del20_ps_outlier.list')
HC_SR = make_ID('SR_HC_ps_outlier.list')
HC_NB = make_ID('NB_HC_ps_outlier.list')

a <- intersect(Del19, Del20)
length(a)
b <- intersect(Del19, HC_SR)
length(b)
c <- intersect(Del19, HC_NB)
length(c)
e <- intersect(Del20, HC_SR)
length(e)
d <- intersect(Del20, HC_NB)
length(d)
f <- intersect(HC_NB, HC_SR)
length(f)

g <- intersect(HC_SR, (intersect(Del19, Del20)))
length(g)

h <- intersect(HC_NB, (intersect(Del19, Del20)))
length(h)

i <- intersect(HC_SR, (intersect(HC_NB, (intersect(Del19, Del20)))))
length(i)

##########################
#  Shared venndiagram    #
##########################
library(ggVennDiagram)
outlier_1<-read.table("Del19_FDR_outlier.list", header=T)
dim(outlier_1)[1]
outlier_2<-read.table("Del20_FDR_outlier.list", header=T)
dim(outlier_2)[1]
outlier_3<-read.table("SR_HC_FDR_outlier.list", header=T)
dim(outlier_3)[1]
outlier_4<-read.table("NB_HC_FDR_outlier.list", header=T)
dim(outlier_4)[1]
x1 = list( Surv_Ref19=outlier_1$id, Surv_Ref_20=outlier_2$id, HC_SR=outlier_3$id, HC_NB=outlier_4$id)

# this dataset has too many SNPs and will cause trouble in venndiagram
# outlier_1<-read.table("Del19_ps_outlier.list", header=T)
# dim(outlier_1)[1]
# outlier_2<-read.table("Del20_ps_outlier.list", header=T)
# dim(outlier_2)[1]
# outlier_3<-read.table("SR_HC_ps_outlier.list", header=T)
# dim(outlier_3)[1]
# outlier_4<-read.table("NB_HC_ps_outlier.list", header=T)
# dim(outlier_4)[1]
# x2 = list(ps_Del19=outlier_1$id, ps_Del20=outlier_2$id, ps_HC_SR=outlier_3$id, ps_HC_NB=outlier_4$id)

jpeg("FDR_venndiagram.jpg", width = 16, height = 9, units = 'in', res = 300)
ggVennDiagram(x1, label_alpha = 0)
dev.off()


#######################
#  Check delta_p      #
#######################

pname = "REF19_CHR19_out_all_z.txt" # 2-side p-value
dat = read.delim(pname, header = FALSE, sep='\t')
# count how many SNPs with positive vs negative frequency change
dat$id = paste0(dat$V1,'_',dat$V2)
outlier <- dat[which(dat$id %in% i),][,1:5]
colnames(outlier) <- c("chr","pos","p1","p2","deltap")
cnt = length(outlier$deltap)
pos_cnt = length(outlier[which(outlier$deltap>0),]$pos)
message("positive change SNPs: ", pos_cnt)
neg_cnt = length(outlier[which(outlier$deltap<0),]$pos)
message("negative change SNPs: ", neg_cnt)
outlier <- outlier[with(outlier, order(-deltap)),]


pname = "Del19_FDR_outlier.list" # 2-side p-value
dat = read.delim(pname, header = TRUE, sep='\t')
length(dat$chr)
outlier=dat[,1:5]
colnames(outlier) <- c("chr","pos","p1","p2","deltap")
cnt = length(outlier$deltap)
pos_cnt = length(outlier[which(outlier$deltap>0),]$pos)
message("positive change SNPs: ", pos_cnt)
neg_cnt = length(outlier[which(outlier$deltap<0),]$pos)
message("negative change SNPs: ", neg_cnt)
outlier <- outlier[with(outlier, order(-deltap)),]

#######################
#  plot delta_p      #
#######################
# plot distribution
jpeg("Comb_Del19_FDR_5870.jpg", width = 16, height = 9, units = 'in', res = 300)
jpeg("Comb_shared_657.jpg", width = 16, height = 9, units = 'in', res = 300)
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

###########################
## test on random shares ##
###########################

n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  o1 = sample(a, 3265)
  o2 = sample(a, 2559)
  cnt = length(intersect(o1,o2))
  boot_cnt[i] = cnt
}

b=t.test(boot_cnt)

hist(boot_cnt)


###########################
## Pearson correlation   ##
###########################
library(ggplot2)
library(grid)
library(gridExtra)

dat1 = read.delim("ps_Del19_challenge.txt", header = FALSE, sep='\t')
length(dat1[which(dat1$V6==0),]$V1)
dat1$V6[dat1$V6 == 0 ] <- 1e-5 
length(dat1[which(dat1$V6==1),]$V1)
dat1$V6[dat1$V6 == 1 ] <- 0.99999
p1 = -log(dat1$V6)
max(p1)
dat2 = read.delim("barnard_REF19_CHR19.txt", header = FALSE, sep='\t') 
p2 = -log(dat2$V6)
max(p2)
df = as.data.frame(cbind(p1, p2))
colnames(df) <- c("SGS_p", "Barnard_p")
head(df)

df <- df[1:100000,]

grob1 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(df$SGS_p, df$Barnard_p), 4) ), x = 0.63, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

jpeg("Del19_pearson.jpg", width = 16, height = 9, units = 'in', res = 300)
ggplot(df, aes(x=SGS_p, y=Barnard_p)) + geom_point() + ggtitle("Del2019 challenge")+
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "SGS -log10(p value)", limits = c(0, 10)) + 
  scale_y_continuous(name = "Barnard's test -log10(p value)", limits = c(0, 10)) + 
  annotation_custom(grob1) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
 
dat1 = read.delim("ps_Del20_challenge.txt", header = FALSE, sep='\t') 
p1 = -log(dat1$V6)
max(p1)
dat2 = read.delim("barnard_REF20_CHR20.txt", header = FALSE, sep='\t') 
p2 = -log(dat2$V6)
max(p2)
df = as.data.frame(cbind(p1, p2))
colnames(df) <- c("SGS_p", "Barnard_p")
head(df)

df <- df[1:1000,]
attach(df)
grob1 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(SGS_p, Barnard_p), 4) ), x = 0.63, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))

jpeg("Del30_pearson.jpg", width = 16, height = 9, units = 'in', res = 300)
ggplot(df, aes(x=SGS_p, y=Barnard_p)) + geom_point() + ggtitle("Del2020 challenge")+
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "SGS -log10(p value)", limits = c(0, 20)) + 
  scale_y_continuous(name = "Barnard's test -log10(p value)", limits = c(0, 20)) + 
  annotation_custom(grob1) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
  