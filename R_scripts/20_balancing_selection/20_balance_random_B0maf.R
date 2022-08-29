#########################################
#######  process outlier for B-score ######
#########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/format/")
###### HC NB outliers ########
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(input, name, distance){
  #input = "ps_Del19_HC_NB.txt"
  dat = read.delim(input, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  message("genome-wide positive ratio is ", length(dat$V5[which(dat$V3-dat$V4>=0)])/length(dat$V5))
  dat$delta_p <- dat$V3-dat$V4
  dat$SNP = paste0(dat$V1,'_',dat$V2)
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
  df1 <- dat[which(dat$adj< 0.05),1:2]
  random = sample(dat$SNP, length(df1$chromo))
  dat_ <- dat[which(dat$SNP %in% random),]
  message(paste0("total number of outliers is ", length(dat_$chromo)))
  message("random snp positive ratio is ", length(dat_$chromo[which(dat_$delta_p>=0)])/length(dat_$chromo))
  dat_ = dat_[with(dat_, order(chromo, position)),]
  bed_list <- paste0(dat_$chromo, "\t" , dat_$position-distance, "\t", dat_$position+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(name, ".random.list"), row.names=F, col.names = F, quote=F, sep="\t")
  return(bed_list)
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random_format/")
for (i in seq(1,100)){
  format_bed("ps_Del19_HC_NB.txt", paste0("./NB_HC_format/NB_HC_random_2507_",i), 250)
}

for (i in seq(1,100)){
  format_bed("ps_Del19_challenge.txt", paste0("./CHR19_REF19_format/CHR19_REF19_random_2932_",i), 250)
}

# challenge_random <- format_bed("ps_Del19_challenge.txt", "CHR19_REF19_2932", 250)
# HC_NB_random <- format_bed("ps_Del19_HC_NB.txt", "NB_HC_2507", 250)
# HC_SR_random <- format_bed("ps_Del19_HC_SR.txt", "SR_HC_2382", 250)

########################################################
# Step 2: using bedtools to merge intervals
########################################################
# path in /workdir/hz269/DelBay_all_angsd_final/15_LD_prunning/process_bed_file

# for f in *.list; do
# echo $f
# cat $f | wc -l
# bedtools merge -i $f > $f'.merged.txt'
# cat $f'.merged.txt' | wc -l
# done

########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################

format_rf <- function(pname){
  dat = read.delim(pname, header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("CHR19_REF19_2932.random.list.merged.txt")
format_rf("NB_HC_2507.random.list.merged.txt")
format_rf("SR_HC_2382.random.list.merged.txt")
# here we will switch to cluster for data running
# the output of angsd run will be used as input for step 5

#######################################
# Step 4: convert mafs to B0maf format 
#######################################

format_mafs <- function(headname){
  dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers.mafs"), header = T)
  #dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers.mafs"), header = T)
  # do formatting for each chromosome
  for (j in c( 'NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j), ]
    # filter out SNPs with maf < 0.05 or maf > 0.95 in each population
    chr_ = c()
    pos_ = c()
    genPos_ = c()
    g_a_ = c()
    g_tot_ = c()
    #for(i in seq(10)){
    for (i in seq(length(DT[, 1]))) {
      chr = DT$chromo[i]
      pos = DT$position[i]
      genPos = "NA"
      sub = "0"
      g_tot = DT$nInd[i] * 2
      g_a1 = round(DT$knownEM[i] * DT$nInd[i] * 2)
      g_a2 = g_tot - g_a1
      if (DT$knownEM[i] < 0.5) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      } else {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a2)
        g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, genPos_, g_a_, g_tot_)
    output = output[order(output$pos_), ]
    print(length(output$pos_))
    write.table(output, file = paste0(headname, "_", j , ".mafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/SR_HC/")
format_mafs("HC")
format_mafs("SR")
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/NB_HC/")
format_mafs("HC")
format_mafs("NB")
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/CHR19_REF19")
format_mafs("REF19")
format_mafs("CHR19")

format_mafs <- function(headname, k){
  dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers_", k, ".mafs"), header = T)
  #dat <- read.table(paste0(headname,"_minmapq30_minq20_CV30_masked_noinvers.mafs"), header = T)
  # do formatting for each chromosome
  for (j in c( 'NC_035780.1', 'NC_035781.1', 'NC_035782.1', 'NC_035783.1', 'NC_035784.1', 'NC_035785.1', 'NC_035786.1', 'NC_035787.1', 'NC_035788.1', 'NC_035789.1')){ 
    DT = dat[which(dat$chromo %in% j), ]
    # filter out SNPs with maf < 0.05 or maf > 0.95 in each population
    chr_ = c()
    pos_ = c()
    genPos_ = c()
    g_a_ = c()
    g_tot_ = c()
    #for(i in seq(10)){
    for (i in seq(length(DT[, 1]))) {
      chr = DT$chromo[i]
      pos = DT$position[i]
      genPos = "NA"
      sub = "0"
      g_tot = DT$nInd[i] * 2
      g_a1 = round(DT$knownEM[i] * DT$nInd[i] * 2)
      g_a2 = g_tot - g_a1
      if (DT$knownEM[i] < 0.5) {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a1)
        g_tot_ = c(g_tot_, g_tot)
      } else {
        pos_ = c(pos_, pos)
        genPos_ = c(genPos_, genPos)
        g_a_ = c(g_a_, g_a2)
        g_tot_ = c(g_tot_, g_tot)
      }
    }
    output = data.frame(pos_, genPos_, g_a_, g_tot_)
    output = output[order(output$pos_), ]               
    print(length(output$pos_))
    write.table(output, file = paste0(headname,"_",k, "_",j, ".mafs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random_format/CHR19_REF19_format")
for (k in seq(1:100)){
  format_mafs("CHR19", k)
} 

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random_format/NB_HC_format/")
for (k in seq(1:100)){
  format_mafs("HC", k)
} 
#####################################
# Step 5: process the output files
#####################################

format_random_output <- function(pop, j, random_list){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  #random_list = "CHR19_REF19_2932.random.list"
  dat_ = read.delim(random_list, header = FALSE, sep='\t')
  dat_$SNP = paste0(dat_$V1,'_',(dat_$V2+dat_$V3)/2)
  
  chr_ = c()
  pos_ = c()
  CLR_ = c()
  nSites_ = c()
  SNP_ = c()
  for (i in seq(length(dat[, 1]))) {
    chr = dat$chromo[i]
    pos = dat$position[i]
    CLR = dat$CLR[i]
    nSites = dat$nSites[i]
    SNP = dat$SNP[i]
    if (dat$SNP[i] %in% dat_$SNP) {
      chr_ = c(chr_, chr)
      pos_ = c(pos_, pos)
      CLR_ = c(CLR_, CLR)
      nSites_ = c(nSites_, nSites)
      SNP_ = c(SNP_, SNP)
    }
  }
  output = data.frame(chr_, pos_, CLR_, nSites_, SNP_)
  colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  output = output[with(output, order(chromo, position)),]
  message(paste0("total number of snp is ", dim(output)[1]))
  #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
  return(output)
}

format_outlier_output <- function(pop, j, outlier_list){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  #outlier_list = "ps_Del19_challenge.txt.outlier.list"
  dat_ = read.delim(outlier_list, header = FALSE, sep='\t')
  dat_$SNP = dat_$V1
  
  chr_ = c()
  pos_ = c()
  CLR_ = c()
  nSites_ = c()
  SNP_ = c()
  for (i in seq(length(dat[, 1]))) {
    chr = dat$chromo[i]
    pos = dat$position[i]
    CLR = dat$CLR[i]
    nSites = dat$nSites[i]
    SNP = dat$SNP[i]
    if (dat$SNP[i] %in% dat_$SNP) {
      chr_ = c(chr_, chr)
      pos_ = c(pos_, pos)
      CLR_ = c(CLR_, CLR)
      nSites_ = c(nSites_, nSites)
      SNP_ = c(SNP_, SNP)
    }
  }
  output = data.frame(chr_, pos_, CLR_, nSites_, SNP_)
  colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  output = output[with(output, order(chromo, position)),]
  message(paste0("total number of snp is ", dim(output)[1]))
  #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
  return(output)
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/random_output/CHR19_REF19")
df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("CHR19", i, "CHR19_REF19_2932.random.list"))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/outlier_output/CHR19_REF19")
df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_outlier_output("CHR19", i, "ps_Del19_challenge.txt.outlier.list"))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]

if(dim(df_1)[1] < dim(df_2)[1]){
  idx = seq(1:dim(df_2)[1])
  random = sort(sample(idx, dim(df_1)[1]))
  df_2 = df_2[random,]
} else {
  idx = seq(1:dim(df_1)[1])
  random = sort(sample(idx, dim(df_2)[1]))
  df_1 = df_1[random,]
}

####### HC_NB random
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/random_output/NB_HC/")
df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("HC", i, "NB_HC_2507.random.list"))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/outlier_output/NB_HC/")
df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_outlier_output("HC", i, "ps_Del19_HC_NB.txt.outlier.list"))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]

if(dim(df_1)[1] < dim(df_2)[1]){
  idx = seq(1:dim(df_2)[1])
  random = sort(sample(idx, dim(df_1)[1]))
  df_2 = df_2[random,]
} else {
  idx = seq(1:dim(df_1)[1])
  random = sort(sample(idx, dim(df_2)[1]))
  df_1 = df_1[random,]
}

####### HC_SR random

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/random_output/SR_HC/")
df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("SR", i, "SR_HC_2382.random.list"))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/outlier_output/SR_HC/")
df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_outlier_output("SR", i, "ps_Del19_HC_SR.txt.outlier.list"))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]

set.seed(200)
if(dim(df_1)[1] < dim(df_2)[1]){
  idx = seq(1:dim(df_2)[1])
  random = sort(sample(idx, dim(df_1)[1]))
  df_2 = df_2[random,]
} else {
  idx = seq(1:dim(df_1)[1])
  random = sort(sample(idx, dim(df_2)[1]))
  df_1 = df_1[random,]
}
#########################
# Step 2: start to plot #
#########################

df_total = rbind(df_1, df_2)
df_total$group = rep(c("outlier", "random"), c(dim(df_1)[1], dim(df_2)[1]))
colnames(df_total)=c('chromo', 'position', 'B0maf_score', 'nSites', 'SNP', 'group')

pwc <- df_total %>% 
  pairwise_t_test(
    B0maf_score ~ group, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

ggboxplot(df_total, x = "group", y = "B0maf_score") +
  #stat_pvalue_manual(pwc, hide.ns = TRUE) +
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  labs( x = NULL, y = expression(Balancing~selection~score~'('~B[0]~','~MAF~')'))+
  theme(text = element_text(size = 15)) 

df2 <- summarySE(df_total, measurevar="B0maf_score", groupvars=c("group"))
df2
ggplot(df2, aes(x=group, y=B0maf_score, colour=group, group=group)) + 
  geom_errorbar(aes(ymin=B0maf_score-se, ymax=B0maf_score+se), colour="black", width=.1, position=pd) +
  coord_cartesian(ylim=c(5.2, 5.8))+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1))+ 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6, hjust = 1.4, vjust = -0)+
  #geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  labs(x = NULL, y = expression(r^2)) +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("OJ", "VC"),
                   labels=c("Orange juice", "Ascorbic acid"),
                   l=40) +                    # Use darker colors, lightness=40
  theme_bw() +
  theme(text = element_text(size = 15)) 

graph2ppt(file="HC_NB.pptx", width=4, height=3)

########################################
######  reformat the data for plot #####
########################################
format_outlier_output <- function(pop, j, outlier_list){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  #outlier_list = "ps_Del19_challenge.txt.outlier.list"
  dat_ = read.delim(outlier_list, header = FALSE, sep='\t')
  dat_$SNP = dat_$V1
  
  chr_ = c()
  pos_ = c()
  CLR_ = c()
  nSites_ = c()
  SNP_ = c()
  for (i in seq(length(dat[, 1]))) {
    chr = dat$chromo[i]
    pos = dat$position[i]
    CLR = dat$CLR[i]
    nSites = dat$nSites[i]
    SNP = dat$SNP[i]
    if (dat$SNP[i] %in% dat_$SNP) {
      chr_ = c(chr_, chr)
      pos_ = c(pos_, pos)
      CLR_ = c(CLR_, CLR)
      nSites_ = c(nSites_, nSites)
      SNP_ = c(SNP_, SNP)
    }
  }
  output = data.frame(chr_, pos_, CLR_, nSites_, SNP_)
  colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  output = output[with(output, order(chromo, position)),]
  message(paste0("total number of snp is ", dim(output)[1]))
  #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
  return(output)
}

format_random_output <- function(pop, j, random_list){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  #random_list = "CHR19_REF19_2932.random.list"
  output = data.frame(dat$chromo, dat$position, dat$CLR, dat$nSites, dat$SNP)
  colnames(output)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
  output = output[with(output, order(chromo, position)),]
  message(paste0("total number of snp is ", dim(output)[1]))
  #write.table(output, paste0(pname, ".outlier.list"), row.names=F, col.names = T, quote=F, sep="\t")
  return(output)
}

####### HC_NB random
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/global_output")
df_1=data.frame()
for (i in seq(0,9)) {
#for (i in seq(4,4)) {
  #df_1 = rbind(df_1,format_outlier_output("CHR19", i, "ps_Del19_challenge.txt.outlier.list"))
  df_1 = rbind(df_1,format_outlier_output("REF19", i, "ps_Del19_challenge.txt.outlier.list"))
  #df_1 = rbind(df_1,format_outlier_output("NB", i, "ps_Del19_HC_NB.txt.outlier.list"))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_2=data.frame()
for (i in seq(0,9)) {
#for (i in seq(4,4)) {
  #df_2 = rbind(df_2,format_random_output("CHR19", i))
  df_2 = rbind(df_2,format_random_output("REF19", i))
  #df_2 = rbind(df_2,format_random_output("NB", i))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

library(dplyr)
library(Rmisc)
library(export)
combine_df <- function(df_1, df_2){
  df_2 = sample_n(df_2, dim(df_1)[1])
  df_1 <- df_1[which(df_1$CLR != 0), ]
  #df_1 <- df_1[which(df_1$nSites >= 1), ]
  df_1 <- df_1[which(df_1$nSites >= 10), ]
  df_2 <- df_2[which(df_2$CLR != 0), ]
  df_2 <- df_2[which(df_2$nSites >= 10), ]
  #df_2 <- df_2[which(df_2$nSites >= 1), ]
  if(dim(df_1)[1] < dim(df_2)[1]){
    idx = seq(1:dim(df_2)[1])
    random = sort(sample(idx, dim(df_1)[1]))
    df_2 = df_2[random,]
  } else {
    idx = seq(1:dim(df_1)[1])
    random = sort(sample(idx, dim(df_2)[1]))
    df_1 = df_1[random,]
  }
  df_total = rbind(df_1, df_2)
  df_total$group = rep(c("outlier", "random"), c(dim(df_1)[1], dim(df_2)[1]))
  return(df_total)
}

random_mean <- c()
for(i in seq(1,1000)){
  df <- combine_df(df_1, df_2)
  df2 <- summarySE(df, measurevar="CLR", groupvars=c("group"))
  message(i)
  print(df2$CLR[2])
  random_mean = c(random_mean, df2$CLR[2])
}

df<-combine_df(df_1, df_2)
observed_mean  <- summarySE(df, measurevar="CLR", groupvars=c("group"))$CLR[1]
random_mean <- as.data.frame(random_mean)
random_mean_point <- mean(random_mean$random_mean)

sum(random_mean$random_mean > observed_mean)

hist_p <- ggplot(random_mean, aes(x=random_mean, fill=NULL)) +
  geom_histogram(fill="white", color="grey", bins = 15)+
  geom_vline(aes(xintercept=mean(random_mean)), color="grey",
             linetype="dashed")+
  labs(x=expression(Balancing~selection~score~"("~B0MAF~")"), y = "Count")+
  theme_classic() +
  ylim(0, 250)+
  theme(text = element_text(size = 15)) +
  geom_text(x=observed_mean-0.0022, label=paste("SGS candidate mean \n", "=", round(observed_mean, digits = 3)), y=15, color="red", size = 4)+
  geom_text(x=random_mean_point+0.0027, label=paste("Mean of null distribution =\n", "=", round(random_mean_point, digits = 3)), y=18, color="Black", size = 4)+
  geom_vline(xintercept=observed_mean, linetype="dashed", color = "red") 
  # xlim(5, 6.5) +
  # ylim(0, 25)

hist_p

graph2ppt(file="HC_NB_NB_B0MAF_comp", width=4, height=3) 
graph2ppt(file="CHR19_REF19_CHR19_B0MAF_comp", width=4, height=3)
graph2ppt(file="CHR19_REF19_REF19_B0MAF_comp", width=4, height=3) 
