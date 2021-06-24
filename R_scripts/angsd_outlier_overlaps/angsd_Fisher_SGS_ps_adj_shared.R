#######################
##  Adjust p-value   ##
#######################
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Indep_Fisher_SGS_shared_FDR_outliers")

# looking at indepedent runs
ps_adj <- function(pname, outname1, outname2) {
  #pname = "fish_REF19_CHR19_REF19_CHR19_NB_HC.txt" # 2-side p-value
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

ps_adj("fish_NB_HC_REF19_CHR19_NB_HC.txt", "fish_NB_HC_REF19_CHR19_NB_HC_ps_outlier.list", "fish_NB_HC_REF19_CHR19_NB_HC_FDR_outlier.list")
ps_adj("fish_REF19_CHR19_REF19_CHR19_NB_HC.txt", "fish_REF19_CHR19_REF19_CHR19_NB_HC_ps_outlier.list", "fish_REF19_CHR19_REF19_CHR19_NB_HC_FDR_outlier.list")

###############################################################
### Bedtools : find the intersection between SNPs and genes ###
###############################################################

# the code is adopted from the tutorial here https://github.com/clairemerot/physalia_adaptation_course/tree/2021/05_day5

###########################################
### format outlier files into bed files ###
###########################################

#######################
#  make bed file      #
#######################
# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Indep_Fisher_SGS_shared_FDR_outliers")

SGS_bed_format <- function(input, output){
  outlier_temp_SGS<-read.table(input, header=T)
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("SGS_Del19_FDR_outlier.list", "./bed/SGS_Del19_FDR_10K.bed" )
bed_format("SGS_NB_HC_FDR_outlier.list", "./bed/SGS_NB_HC_FDR_10K.bed" )


###########################################
#  make bed file for Fisher outliers      #
###########################################

# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Indep_Fisher_SGS_shared_FDR_outliers")

fisher_bed_format <- function(input, output){
  outlier_temp_SGS<-read.table(input, header=T)
  #outlier_temp_SGS<-read.table("Fisher_NB_HC_REF19_CHR19_NB_HC_FDR_outlier.list", header=T)
  #what's its dimension?
  dim(outlier_temp_SGS)
  #which size around the SNP
  #window<-2000
  window<-10000
  #add a vector with start position
  outlier_temp_SGS$start<-outlier_temp_SGS$pos-(window/2)
  #oups it can't be ngative! replace negative by 0
  outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
  #add a vector with stop position
  outlier_temp_SGS$stop<-outlier_temp_SGS$pos+(window/2)
  #have a look
  head(outlier_temp_SGS)
  #which columns shoud we keep?
  outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,9,10,8)]
  colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
  head(outlier_temp_SGS_bed)
  #save your file
  write.table(outlier_temp_SGS_bed, output, row.names=F, sep="\t", quote=F,col.names=F)
}

bed_format("Fisher_NB_HC_FDR_outlier.txt", "./bed/Fisher_NB_HC_10K.bed" )
bed_format("Fisher_Del19_FDR_outlier.txt", "./bed/Fisher_Del19_10K.bed" )

###########################################
#       merge the intersect regions       #
###########################################
#====================MERGE======================
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Indep_Fisher_SGS_shared_FDR_outliers/bed")

merge <- function(input, output){
  #outlier_DT = read.delim("Fisher_Del19_10K.bed", header = FALSE, sep='\t')
  outlier_DT = read.delim(input, header = FALSE, sep='\t')
  colnames(outlier_DT) <- c("Chromosome", "Start", "End", "id")
  outlier_DT = outlier_DT[order(outlier_DT$Chromosome, outlier_DT$Start),]
  library(hash)
  dict = hash()
  for(i in seq(length(outlier_DT[,1]))){
    chr = outlier_DT$Chromosome[i]
    st = outlier_DT$Start[i]
    ed = outlier_DT$End[i]
    ranges = dict[[chr]]
    if(length(ranges)>0 && ranges[[length(ranges)]][2]>=st){ #overlap
      last_range = ranges[[length(ranges)]]
      ranges = ranges[-length(ranges)] #remove last one
      dict[[chr]] = append(ranges, list(c(last_range[1], max(last_range[2], ed))))
    }else{
      dict[[chr]] = append(ranges, list(c(st, ed)))
      #previous zFst set calculate median, write,. reset zFst set
    }
  }
  #===================================================
  CHR = c()
  ST = c()
  ED = c()
  for(chr in unique(outlier_DT$Chromosome)){
    ranges = dict[[chr]] 
    for(i in seq(length(ranges))){
      CHR = c(CHR, chr)
      ST = c(ST, ranges[[i]][1])
      ED = c(ED, ranges[[i]][2])
    }
  }
  outlier_DT = data.frame(Chromosome=CHR, Start=ST, End=ED, Feature="merged")
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature')
  write.table(outlier_DT, file = paste0(output, ".merged.igv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

merge("Fisher_Del19_10K.bed", "Fisher_Del19_10K.bed")
merge("Fisher_NB_HC_10K.bed", "Fisher_NB_HC_10K.bed")
merge("SGS_Del19_FDR_10K.bed", "SGS_Del19_FDR_10K.bed")
merge("SGS_NB_HC_FDR_10K.bed", "SGS_NB_HC_FDR_10K.bed")


#merge("SGS_Fisher_NB_HC_FDR_10K.intersect", "SGS_Fisher_NB_HC_FDR_10K.intersect")
#merge("SGS_Fisher_REF19_CHR19_FDR_10K.intersect", "SGS_Fisher_REF19_CHR19_FDR_10K.intersect")

###########################################
#       merge the intersect regions       #
###########################################

setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Fisher_SGS_shared_FDR_outliers/bed")
system("bedtools intersect -a Fisher_Del19_10K.bed.merged.igv -b SGS_Del19_FDR_10K.bed.merged.igv > SGS_Fisher_REF19_CHR19_FDR_10K.intersect", intern = TRUE)
system("bedtools intersect -a Fisher_NB_HC_10K.bed.merged.igv -b SGS_NB_HC_FDR_10K.bed.merged.igv > SGS_Fisher_NB_HC_FDR_10K.intersect", intern = TRUE)
system("bedtools intersect -a SGS_Fisher_REF19_CHR19_FDR_10K.intersect -b SGS_Fisher_NB_HC_FDR_10K.intersect > SGS_Fisher_shared_10K_merged.intersect", intern = TRUE)

####################################################################################################################################
#######################
#  Shared outliers    #
#######################
setwd("/Volumes/cornell/Cohort_adaptation/DelBay_adult/Fisher_SGS_shared_FDR_outliers")
make_ID <- function(file_name){
  #file_name ="Del19_outlier.list"
  dat = read.delim(file_name, header = TRUE, sep='\t')
  dat = dat[with(dat, order(chr, pos)),]
  message(length(dat$id))
  return(dat$id)
}

SGS_Del19 = make_ID('SGS_Del19_FDR_outlier.list')
Fisher_Del19 = make_ID("Fisher_REF19_CHR19_REF19_CHR19_NB_HC_FDR_outlier.list")


Del20 = make_ID('Del20_FDR_outlier.list')
HC_SR = make_ID('SR_HC_FDR_outlier.list')
HC_NB = make_ID('NB_HC_FDR_outlier.list')

a <- intersect(SGS_Del19, Fisher_Del19)
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
  