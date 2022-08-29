########################################
# Measure the SNP density for outliers #
########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/window_snp_count/")
dat = read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", header = FALSE, sep='\t')
snp_all <- data.frame(dat$V1, dat$V2, dat$V2)
write.table(snp_all, "snp_all.bed", row.names=F, col.names = F, quote=F, sep="\t")

bedtools = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/bedtools2/bin/bedtools"
system(paste(bedtools," intersect -a 18_HC_NB_shared_outlier.bed -b snp_all.bed -c > 18_HC_NB.count.txt", sep=""))
system(paste(bedtools," intersect -a 19_HC_NB_shared_outlier.bed -b snp_all.bed -c > 19_HC_NB.count.txt", sep=""))
system(paste(bedtools," intersect -a 21_HC_NB_shared_outlier.bed -b snp_all.bed -c > 21_HC_NB.count.txt", sep=""))
system(paste(bedtools," intersect -a 19_Sur_Ref_shared_outlier.bed -b snp_all.bed -c > 19_Sur_Ref.count.txt", sep=""))
system(paste(bedtools," intersect -a 20_Sur_Ref_shared_outlier.bed -b snp_all.bed -c > 20_Sur_Ref.count.txt", sep=""))
system(paste(bedtools," intersect -a CD5_outlier.bed.merged.txt -b snp_all.bed -c > CD5_outlier.bed.merged.count.txt", sep=""))
system(paste(bedtools," intersect -a 18_HC_NB_random_2295_1.bed.merged.txt -b snp_all.bed -c > 18_HC_NB_random_2295_1.bed.merged.count.txt", sep=""))
system(paste(bedtools," intersect -a CD5_random_119_2.bed.merged.txt -b snp_all.bed -c > CD5_random_119_2.bed.merged.count.txt", sep=""))
system(paste(bedtools," intersect -a 19_HC_NB_random_8622_2.bed.merged.txt -b snp_all.bed -c > 19_HC_NB_random_8622_2.bed.merged.count.txt", sep=""))



data = read.delim("18_HC_NB.count.txt", header = FALSE, sep='\t')
mean1 <- mean(data$V4)
data = read.delim("19_HC_NB.count.txt", header = FALSE, sep='\t')
mean2 <- mean(data$V4)
data = read.delim("21_HC_NB.count.txt", header = FALSE, sep='\t')
mean3 <- mean(data$V4)
data = read.delim("19_Sur_Ref.count.txt", header = FALSE, sep='\t')
mean4 <- mean(data$V4)
data = read.delim("20_Sur_Ref.count.txt", header = FALSE, sep='\t')
mean5 <- mean(data$V4)
data = read.delim("CD5_outlier.bed.merged.count.txt", header = FALSE, sep='\t')
mean6 <- mean(data$V4)
data = read.delim('18_HC_NB_random_2295_1.bed.merged.count.txt', header = FALSE, sep='\t')
mean7 <- mean(data$V4)


# on average how many snps within a 500 bp window?
mean(mean1, mean2, mean3, mean4, mean5)
# 12


###################################################################
# format the total outliers and exclude them from random sampling #
###################################################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/shared_outliers")
# cat *.txt > total_outliers.txt
dat = read.delim("total_outliers.txt", header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
write.table(dat, "total_outliers.sort.txt", row.names=F, col.names = F, quote=F, sep="\t")
#./merge.sh
# find the SNPs within outlier intervals
dat1 <- read.delim("total_outliers.sort.txt.merged.txt", header = FALSE, sep='\t')
dat1_ <- paste0(dat1$V1, "\t" , dat1$V2-500, "\t", dat1$V3+500)
write.table(dat1_, "total_outlier.sort.merged.buffer.bed", row.names=F, col.names = F, quote=F, sep="\t")

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/window_snp_count/")
dat = read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", header = FALSE, sep='\t')
snp_all <- data.frame(dat$V1, dat$V2, dat$V2)
write.table(snp_all, "snp_all.bed", row.names=F, col.names = F, quote=F, sep="\t")

bedtools = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/bedtools2/bin/bedtools"
system(paste(bedtools," intersect -a total_outlier.sort.merged.buffer.bed -b snp_all.bed -c > total_outlier.sort.merged.buffer.count.txt", sep=""))
system(paste(bedtools, " subtract -a snp_all.bed -b total_outlier.sort.merged.buffer.bed > neutral_snp_list.txt", sep=""))

#############################################################
#########  output the outlier and random SNP id   ##########
#############################################################
# 
# colnames(snp_all) = c("chromo", "position", "position")
# dat = snp_all
# # this block of code is used to elimiate SNPs with flanking sites less than a certain amount of numbers (for example, within a 500 bp window a focal SNP must have > 10 SNPs in the window)
# idx = dat$position>0 # all TRUE
# win = 250
# 
# #dat= dat[1:10000,]
# #idx = idx[1:10000]
# 
# for(i in seq(dim(dat)[1])){
# #for(i in seq(1:10000)){
#   chr = dat$chromo[i]
#   pos = dat$position[i]
#   cnt = 0
#   # pre range
#   for(j in seq(i,i-20)){
#     if(j<1||j>dim(dat)[1])
#       break
#     cur_chr = dat$chromo[j]
#     cur_pos = dat$position[j]
#     if(cur_chr==chr && cur_pos>(pos-win))
#       cnt = cnt + 1
#   }
#   # post range
#   for(j in seq(i,i+20)){
#     if(j<1||j>dim(dat)[1])
#       break
#     cur_chr = dat$chromo[j]
#     cur_pos = dat$position[j]
#     if(cur_chr==chr && cur_pos<(pos+win))
#       cnt = cnt + 1
#   }
#   if(cnt<12)
#     idx[i] = FALSE
# }
# 
# dat_min12 = dat[idx,]

# n = 1178698 remainted for random test
#df1 <- dat[which(dat$adj< 0.05),1:2]
#df2 <- paste0(df1$chromo, "\t", df1$position, "\t")
#message(paste0("number of outlier in ", pname, " is ", length(df1$chromo)))
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
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD")
dat = read.delim("neutral_snp_list.txt", header = FALSE, sep='\t')
colnames(dat) = c("chromo", "position", "position")

format_bed <- function(dat, cnt, name, distance){
  dat$SNP <- paste0(dat$chromo, dat$position)
  random = sample(dat$SNP, cnt)
  df_ <- dat[which(dat$SNP %in% random),1:2]
  df_ = df_[with(df_, order(chromo, position)),]
  bed_list <- paste0(df_$chromo, "\t" , df_$position-distance, "\t", df_$position+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(name, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

for (i in seq(1,100)){
  format_bed(dat, 2295, paste0("./18_HC_NB/18_HC_NB_random_2295_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 2147, paste0("./19_HC_NB/19_HC_NB_random_2147_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 1893, paste0("./21_HC_NB/21_HC_NB_random_1893_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 2250, paste0("./19_Sur_Ref/19_Sur_Ref_random_2250_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 2294, paste0("./20_Sur_Ref/20_Sur_Ref_random_2294_",i), 250)
}

for (i in seq(1,100)){
  format_bed(dat, 119, paste0("./CD5/CD5_random_119_",i), 250)
}

########################################################
# Step 2: using bedtools to merge intervals
########################################################
# merge.sh in /Users/HG/Dropbox/Mac/Documents/HG/DelBay_final/15_local_LD/20_Sur_Ref

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

for (i in seq(1,100)){
  format_rf(paste0("./18_HC_NB/18_HC_NB_random_2295_", i, ".bed.merged.txt"), 1966)
}

for (i in seq(1,100)){
  format_rf(paste0("./19_HC_NB/19_HC_NB_random_2147_", i, ".bed.merged.txt"), 1877)
}

for (i in seq(1,100)){
  format_rf(paste0("./21_HC_NB/21_HC_NB_random_1893_", i, ".bed.merged.txt"), 1583)
}

for (i in seq(1,100)){
  format_rf(paste0("./19_Sur_Ref/19_Sur_Ref_random_2250_", i, ".bed.merged.txt"), 1957)
}

for (i in seq(1,100)){
  format_rf(paste0("./20_Sur_Ref/20_Sur_Ref_random_2294_", i, ".bed.merged.txt"), 1925)
}

for (i in seq(1,100)){
  format_rf(paste0("./CD5/CD5_random_119_", i, ".bed.merged.txt"), 115)
}

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
