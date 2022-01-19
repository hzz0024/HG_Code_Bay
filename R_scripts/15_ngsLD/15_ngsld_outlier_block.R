#############################################
## reorder the ps file from SGS test       ## 
############################################# 

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS")
########### be really careful about the issue from order. clear everthing in the data before doing that.
pname = "ps_CHR19_REF19.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
dat = dat[order(dat[,1], dat[,2]),]
write.table(dat, "./ps_CHR19_REF19_reorder.txt", row.names=F, col.names=F, quote=F, sep="\t")

####################################
##########  output snp list ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay_all_angsd_final/11_SGS/no_shared/indep_SGS_ps")
# output global snp list
df = data.frame(dat$V1, dat$V2)
write.table(df, "./SNP_list/SGS_global_snp_list.txt", row.names=F, col.names = F, quote=F, sep="\t")
# output outliers with FDR < 0.05
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'delta_p', 'ps', 'raw_candidates', 'adj')
dat1 <- dat[which(dat$adj< 0.05),]
dat1$id <- paste0(dat1$chromo,'_',dat1$position)
dim(dat1)
dat1 = dat1[with(dat1, order(chromo, position)),]
df1 <- data.frame(dat1$chromo, dat1$position)
write.table(df1, "./SNP_list/CHR19_REF19_SGS_FDR05.txt", row.names=F, col.names = F, quote=F, sep="\t")

# dat1[which(dat1$position == '16552716'),] check if this chr5_16552716 still there, answer is yes!

#########################################
#######  process outlier for ngsLD ######
#########################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/15_ngsLD/outliter_block_format/")
########################################################
# Step 1: format the ps file and convert it to bed format
########################################################
# must sort the position for ngsLD running
format_bed <- function(pname, distance){
  pname = "ps_CHR20_REF20.txt"
  dat = read.delim(pname, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  message("genome-wide positive ratio is ", length(dat$V5[which(dat$V3-dat$V4>=0)])/length(dat$V5))
  dat$delta_p <- dat$V3-dat$V4
  dat$V6[dat$V6 == 0] = 0.00001
  dat$adj = p.adjust(dat$V6, method = 'BH')
  colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'at_D', 'ps', 'raw_candidates', 'delta_p', 'adj')
  dat_ <- dat[which(dat$adj< 0.05),] # FDR < 0.05
  message(paste0("total number of outliers is ", length(dat_$chromo)))
  message("outlier positive ratio is ", length(dat_$chromo[which(dat_$delta_p>=0)])/length(dat_$chromo))
  dat_ = dat_[with(dat_, order(chromo, position)),]
  bed_list <- paste0(dat_$chromo, "\t" , dat_$position-distance, "\t", dat_$position+distance)
  print(head(bed_list))
  write.table(bed_list, paste0(pname, ".bed"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_bed("ps_CHR20_REF20.txt", 250)
format_bed("ps_Del19_challenge.txt", 250)
format_bed("ps_Del19_HC_NB.txt", 250)
format_bed("ps_Del19_HC_SR.txt", 250)
format_bed("ps_Del19_ARN_COH.txt", 250)
format_bed("ps_Del19_REF19_SR.txt", 250)
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

# ps_Del19_ARN_COH.txt.bed
# 1113
# 922
# ps_Del19_challenge.txt.bed
# 2185
# 1854
# ps_Del19_HC_NB.txt.bed
# 2073
# 1696
# ps_Del20_challenge.txt.bed
# 3117
# 2525

########################################################
# Step 3: convert the bed format to rf input for Angsd
########################################################

format_rf <- function(pname){
  dat = read.delim(pname, header = FALSE, sep='\t')
  angsd_list <- paste0(dat$V1, ":" , dat$V2, "-", dat$V3)
  write.table(angsd_list, paste0(strsplit(pname, split = ".txt")[[1]][1], ".rf.txt"), row.names=F, col.names = F, quote=F, sep="\t")
}

format_rf("ps_Del19_ARN_COH.txt.bed.merged.txt")
format_rf("ps_Del19_REF19_SR.txt.bed.merged.txt")

format_rf("ps_Del19_challenge.txt.bed.merged.txt")
format_rf("ps_Del19_HC_NB.txt.bed.merged.txt")
format_rf("ps_Del19_HC_SR.txt.bed.merged.txt")

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
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_HC_NB.rf.txt > 'ps_Del19_HC_NB.rf.chr'$j'.txt'
# done
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_challenge.rf.txt > 'ps_Del19_challenge.rf.chr'$j'.txt'
# done
# 
# for i in {0..9}; do
# #for i in NC_035780.1 NC_035781.1 NC_035782.1 NC_035783.1 NC_035784.1 NC_035785.1 NC_035786.1 NC_035787.1 NC_035788.1 NC_035789.1; do
# j=$[i+1]
# grep 'NC_03578'$i'.1' ps_Del19_ARN_COH.rf.txt > 'ps_Del19_ARN_COH.rf.chr'$j'.txt'
# done
