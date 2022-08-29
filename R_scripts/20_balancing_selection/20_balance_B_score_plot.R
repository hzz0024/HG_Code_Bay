#########################################
#######  process outlier for B-score ######
#########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/outlier_output/")
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname){
  #pname = "ps_Del19_HC_NB.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  message("genome-wide positive ratio is ", length(dat$V5[which(dat$V3-dat$V4>=0)])/length(dat$V5))
  dat$delta_p <- dat$V3-dat$V4
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
  dat_ <- dat[which(dat$adj< 0.05),] # FDR < 0.05
  message(paste0("total number of outliers is ", length(dat_$chromo)))
  message("outlier positive ratio is ", length(dat_$chromo[which(dat_$delta_p>=0)])/length(dat_$chromo))
  dat_ = dat_[with(dat_, order(chromo, position)),]
  bed_list <- paste0(dat_$chromo, "_" , dat_$position)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".outlier.list"), row.names=F, col.names = F, quote=F, sep="\t")
  return(bed_list)
}

challenge_outliter <- format_bed("ps_Del19_challenge.txt")
HC_NB_outlier <- format_bed("ps_Del19_HC_NB.txt")
HC_SR_outlier <- format_bed("ps_Del19_HC_SR.txt")
HC_SR_outlier <- format_bed("")


format_output <- function(pop, j){
    dat <- read.table(paste0("./CHR19_REF19/",pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
    dat$chromo = paste0("NC_03578", j, '.1')
    dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
    dat = dat[with(dat, order(chromo, physPos)),]
    colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
    
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
      if (dat$SNP[i] %in% challenge_outliter) {
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


df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_output("CHR19", i))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}


df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("REF19", i))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]
df_1 = df_1[which(df_1$SNP %in% df_2$SNP),]
df_2 = df_2[which(df_2$SNP %in% df_1$SNP),]
dim(df_1)
#########################
# Step 2: start to plot #
#########################

df_total = rbind(df_1, df_2)
df_total$pop = rep(c("CHR19", "REF19"), c(dim(df_1)[1], dim(df_2)[1]))
colnames(df_total)=c('chromo', 'position', 'B0maf_score', 'nSites', 'SNP', 'pop')

p <- ggplot(df_total, aes(x=pop, y=B0maf_score, fill=pop)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw()+
  theme(text = element_text(size=30)) 
p

t.test(df_1$CLR, df_2$CLR)






###### HC NB outliers ########

########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname){
  #pname = "ps_Del19_HC_NB.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  message("genome-wide positive ratio is ", length(dat$V5[which(dat$V3-dat$V4>=0)])/length(dat$V5))
  dat$delta_p <- dat$V3-dat$V4
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
  dat_ <- dat[which(dat$adj< 0.05),] # FDR < 0.05
  message(paste0("total number of outliers is ", length(dat_$chromo)))
  message("outlier positive ratio is ", length(dat_$chromo[which(dat_$delta_p>=0)])/length(dat_$chromo))
  dat_ = dat_[with(dat_, order(chromo, position)),]
  bed_list <- paste0(dat_$chromo, "_" , dat_$position)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".outlier.list"), row.names=F, col.names = F, quote=F, sep="\t")
  return(bed_list)
}

challenge_outliter <- format_bed("ps_Del19_challenge.txt")
HC_NB_outlier <- format_bed("ps_Del19_HC_NB.txt")
HC_SR_outlier <- format_bed("ps_Del19_HC_SR.txt")



format_output <- function(pop, j){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  
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
    if (dat$SNP[i] %in% HC_NB_outlier) {
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


df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_output("HC", i))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}


df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("NB", i))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]
df_1 = df_1[which(df_1$SNP %in% df_2$SNP),]
df_2 = df_2[which(df_2$SNP %in% df_1$SNP),]
dim(df_1)
#########################
# Step 2: start to plot #
#########################

df_total = rbind(df_1, df_2)
df_total$pop = rep(c("HC", "NB"), c(dim(df_1)[1], dim(df_2)[1]))
colnames(df_total)=c('chromo', 'position', 'B0maf_score', 'nSites', 'SNP', 'pop')

p <- ggplot(df_total, aes(x=pop, y=B0maf_score, fill=pop)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw()+
  theme(text = element_text(size=30)) 
p

t.test(df_1$CLR, df_2$CLR)





###############  HC_SR_outlier ###############


########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname){
  #pname = "ps_Del19_HC_NB.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  message("genome-wide positive ratio is ", length(dat$V5[which(dat$V3-dat$V4>=0)])/length(dat$V5))
  dat$delta_p <- dat$V3-dat$V4
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
  dat_ <- dat[which(dat$adj< 0.05),] # FDR < 0.05
  message(paste0("total number of outliers is ", length(dat_$chromo)))
  message("outlier positive ratio is ", length(dat_$chromo[which(dat_$delta_p>=0)])/length(dat_$chromo))
  dat_ = dat_[with(dat_, order(chromo, position)),]
  bed_list <- paste0(dat_$chromo, "_" , dat_$position)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".outlier.list"), row.names=F, col.names = F, quote=F, sep="\t")
  return(bed_list)
}

challenge_outliter <- format_bed("ps_Del19_challenge.txt")
HC_NB_outlier <- format_bed("ps_Del19_HC_NB.txt")
HC_SR_outlier <- format_bed("ps_Del19_HC_SR.txt")



format_output <- function(pop, j){
  dat <- read.table(paste0(pop, "_NC_03578",j,".1.mafs.output.500bp.s1.txt"), header = T)
  dat$chromo = paste0("NC_03578", j, '.1')
  dat$SNP <- paste0(dat$chromo,'_',dat$physPos)
  dat = dat[with(dat, order(chromo, physPos)),]
  colnames(dat)=c('position', 'genPos', 'CLR', 'x_hat', 's_hat', 'A_hat',  'nSites', 'chromo', 'SNP')
  
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
    if (dat$SNP[i] %in% HC_SR_outlier) {
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


df_1=data.frame()
for (i in seq(0,9)) {
  df_1 = rbind(df_1,format_output("HC", i))
  colnames(df_1)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_2=data.frame()
for (i in seq(0,9)) {
  df_2 = rbind(df_2,format_output("SR", i))
  colnames(df_2)=c('chromo', 'position', 'CLR', 'nSites', 'SNP')
}

df_1 <- df_1[which(df_1$CLR != 0), ]
df_1 <- df_1[which(df_1$nSites >= 10), ]
df_2 <- df_2[which(df_2$CLR != 0), ]
df_2 <- df_2[which(df_2$nSites >= 10), ]
df_1 = df_1[which(df_1$SNP %in% df_2$SNP),]
df_2 = df_2[which(df_2$SNP %in% df_1$SNP),]
dim(df_1)
#########################
# Step 2: start to plot #
#########################

df_total = rbind(df_1, df_2)
df_total$pop = rep(c("HC", "SR"), c(dim(df_1)[1], dim(df_2)[1]))
colnames(df_total)=c('chromo', 'position', 'B0maf_score', 'nSites', 'SNP', 'pop')

p <- ggplot(df_total, aes(x=pop, y=B0maf_score, fill=pop)) +
  geom_boxplot(alpha=0.7) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1)) + 
  stat_summary(fun = mean, color = "black", geom ="point", aes(group = 1), size = 2, show.legend = FALSE)+
  stat_summary(aes(label=round(..y..,4)), fun=mean, geom="text", size=6,vjust = -0.5)+
  guides(fill = "none") +
  theme_bw()+
  theme(text = element_text(size=30)) 
p

t.test(df_1$CLR, df_2$CLR)
