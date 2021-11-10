#############################################################
#########  output the outlier and random SNP id   ##########
#############################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/random_SNPs_nomatch")
pname = "ps_Del19_HC_NB.txt"
pname = "ps_Del19_challenge.txt"
pname = "ps_Del19_ARN_COH.txt"
pname = "ps_Del19_REF19_SR.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
# process outliers
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
head(dat)

# this block of code is used to elimiate SNPs with flanking sites less than a certain amount of numbers (for example, within a 500 bp window a focal SNP must have > 10 SNPs in the window)
idx = dat$position>0 # all TRUE
win = 250
for(i in seq(dim(dat)[1])){
  #for(i in seq(1:10000)){
  chr = dat$chromo[i]
  pos = dat$position[i]
  cnt = 0
  # pre range
  for(j in seq(i,i-10)){
    if(j<1||j>dim(dat)[1])
      break
    cur_chr = dat$chromo[j]
    cur_pos = dat$position[j]
    if(cur_chr==chr && cur_pos>(pos-win))
      cnt = cnt + 1
  }
  # post range
  for(j in seq(i,i+10)){
    if(j<1||j>dim(dat)[1])
      break
    cur_chr = dat$chromo[j]
    cur_pos = dat$position[j]
    if(cur_chr==chr && cur_pos<(pos+win))
      cnt = cnt + 1
  }
  if(cnt<13)
    idx[i] = FALSE
}

dat1 =dat[idx,]


df1 <- dat[which(dat$adj< 0.05),1:2]
df2 <- paste0(df1$chromo, "\t", df1$position, "\t")
message(paste0("number of outlier in ", pname, " is ", length(df1$chromo)))
# random sample snps from the genome
# random = sample(dat$SNP, length(df1$chromo))
# df3 <- dat[which(dat$SNP %in% random),1:2]
#write.table(df3, "./ps_Del19_HC_NB_random.txt", row.names=F, col.names=F, quote=F, sep="\t")

# check the distribution of p0 p1
#hist(dat[which(dat$adj< 0.05),]$p0)
#hist(dat[which(dat$SNP %in% random),]$p0)

################################################################
### step 1 format the random snps and convert to bed file ######
################################################################

format_bed <- function(dat, cnt, name, distance){
  random = sample(dat$SNP, cnt)
  df3 <- dat[which(dat$SNP %in% random),1:2]
  df3 = df3[with(df3, order(chromo, position)),]
  bed_list <- paste0(df3$chromo, "\t" , df3$position-distance, "\t", df3$position+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(name, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

#format_bed(df3, "REF19_CHR19_random_1598", 250)
#format_bed(df3, "NB_HC_random_1103", 250)

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/zoom_in")
for (i in seq(1,100)){
  format_bed(dat1, 1, paste0("./format/REF19_CHR19_random_1_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 2507, paste0("./NB_HC_format/NB_HC_random_2507_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 2932, paste0("./CHR19_REF19_format/REF19_CHR19_random_2932_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 1507, paste0("./ARN_COH_format/ARN_COH_random_1507_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 999, paste0("./REF19_SR_format/REF19_SR_random_999_",i), 250)
}



########################################################
# Step 2: using bedtools to merge intervals
########################################################
# path in /workdir/hz269/DelBay_all_angsd_final/15_LD_prunning/process_bed_file

# for f in *.bed; do
# echo $f
# cat $f | wc -l
# bedtools merge -i $f > $f'.merged.txt'
# cat $f'.merged.txt' | wc -l
# done

########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################

format_rf <- function(pname, block_number){
  #pname = "SR_HC_format/SR_HC_random_2382.bed.merged.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat$angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  random = sample(dat$angsd_list, block_number)
  dat = dat[which(dat$angsd_list %in% random),]
  avg_leg <- mean(dat$V3-dat$V2)
  message(paste0("average length of block is ", avg_leg))
  write.table(dat$angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("NB_HC_random_1103.bed.merged.txt", 951)
format_rf("REF19_CHR19_random_1598.bed.merged.txt", 1377)
format_rf("SR_HC_random_1051.bed.merged.txt", 926)

format_rf("NB_HC_format/NB_HC_random_2507.bed.merged.txt", 2109)
format_rf("CHR19_REF19_format/CHR19_REF19_random_2932.bed.merged.txt", 2498)
format_rf("SR_HC_format/SR_HC_random_2382.bed.merged.txt", 2025)

########################################################
# Step 4: extract the rf part based on chromosomes
########################################################
# path /Users/HG/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/shared/ngsld_outlier_block
# path /Users/HG/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/shared/random_theta_distribution
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del20_challenge.rf.txt > 'ps_Del20_challenge.rf.chr'$j'.txt'
# done

################################################################
#########  plot the theta distribution for outliers   ##########
################################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/ngsld_theta_distribution")
pname = "test.output.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
colnames(dat) = c('win', 'Chr', 'WinCenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'Tajima', 'fuf', 'fud', 'fayh', 'zeng', 'nSites')
dat$theta <- dat$tP/dat$nSites
hist(dat$theta, probability = T)
lines(density(dat$theta), col=2)
# test for Normality
ks.test(dat$theta, "pnorm", mean=mean(dat$theta), sd=sd(dat$theta))
# find the 95% CI
summary(dat$theta)
t.test(dat$theta, conf.level = 0.95)

mean(dat$theta)
sd(dat$theta)

#########################################################################
#########  random sample snps from the genome with similar pi  ##########
#########################################################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/no_shared/theta_distribution")
pname = "NB.all.snps.theta.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
colnames(dat) = c('win', 'Chr', 'WinCenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'Tajima', 'fuf', 'fud', 'fayh', 'zeng', 'nSites')
dat$theta <- dat$tP/dat$nSites
dat$id <- paste0(dat$Chr, '_', dat$WinCenter)
dat <- dat[complete.cases(dat), ]
# exclude the outliers
dat <- subset(dat, !(dat$id %in% outlier))
hist(dat$theta)
# test for Normality
ks.test(dat$theta, "pnorm", mean=mean(dat$theta), sd=sd(dat$theta))
# find the 95% CI
summary(dat$theta)
t.test(dat$theta, conf.level = 0.95)

# total amount of snps with theta in the target interval
dim(dat[which(dat$theta > a$conf.int[1] & dat$theta < a$conf.int[2]),])

res_dat <- dat[0,]
idx=seq(dim(dat)[1])
while(nrow(res_dat)<2000){
  id=sample(idx,1)
  idx=idx[idx!=id]
  val = dat$theta[id]
  if(val > a$conf.int[1] & val < a$conf.int[2]){
    print(val)
    res_dat = rbind(res_dat, dat[id,])
  }	
}

hist(res_dat$theta)
mean(res_dat$theta)
res_dat$id = paste0(res_dat$Chr,'_',res_dat$WinCenter)

# double check if the random snps are in the global list
# setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/15_ngsLD/ngsld_theta_distribution/NB")
# # extrac the snp id from origianl ps file
# pname = "ps_Del19_HC_NB.txt"
# ps_dat = read.delim(pname, header = FALSE, sep='\t')
# ps_dat = ps_dat[with(ps_dat, order(V1, V2)),]
# ps_dat$id <- paste0(ps_dat$V1, '_', ps_dat$V2)
# colnames(ps_dat)=c('chromo', 'position', 'p1', 'p0', 'dp', 'ps', 'raw_candidates', 'id')
# ps_random = ps_dat[which(res_dat$id %in% ps_dat$id),]

#########################################################
### format the random snps and convert to bed file ######
#########################################################

format_bed <- function(res_dat, name, distance){
  res_dat = res_dat[order(res_dat$Chr, res_dat$WinCenter),]
  bed_list <- paste0(res_dat$Chr, "\t" , res_dat$WinCenter-distance, "\t", res_dat$WinCenter+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(name, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed(res_dat, "NB_random_2000", 500)
