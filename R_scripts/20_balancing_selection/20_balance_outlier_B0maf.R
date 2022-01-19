#########################################
#######  process outlier for B-score ######
#########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/B_score/outlier_output/")
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname, distance){
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
  bed_list_ <- paste0(dat_$chromo, "\t" , dat_$position-distance, "\t", dat_$position+distance)
  print(head(bed_list_))
  write.table(bed_list_, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
  return(bed_list)
}

challenge_outliter <- format_bed("ps_Del19_challenge.txt", 250)
HC_NB_outlier <- format_bed("ps_Del19_HC_NB.txt", 250)
HC_SR_outlier <- format_bed("ps_Del19_HC_SR.txt", 250)
HC_COH_outlier <- format_bed("ps_Del19_HC_COH.txt", 250)


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

format_rf("ps_Del19_HC_COH.txt.bed.merged.txt")
# here we will switch to cluster for data running
# the output of angsd run will be used as input for step 5

#######################################
# Step 4: convert mafs to B0maf format 
#######################################

format_mafs <- function(headname){
  dat <- read.table(paste0(headname,"_all_minq20_minmq30_CV30_masked.mafs"), header = T)
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

setwd("~/Dropbox/Mac/Documents/HG/DelBay20_adult/20_balancing_selection/format")
format_mafs("REF20")
format_mafs("CHR20")

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/SR_HC/")
format_mafs("HC")
format_mafs("SR")
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/NB_HC/")
format_mafs("HC")
format_mafs("NB")
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/20_balancing_selection/random/CHR19_REF19")
format_mafs("REF19")
format_mafs("CHR19")

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
