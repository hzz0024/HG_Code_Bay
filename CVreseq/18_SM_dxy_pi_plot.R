# python script used for dxy estimating is https://github.com/hugang123/Dxy
# for i in {1..10}; do
#   for pop in SL_OBOYS2 CS_NEH CS_DEBY; do
#   vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.vcf --chr $i --keep  $pop --recode --recode-INFO-all --out $pop'_chr'$i
#   done
# done
# usage: win_Dxy.py [-h] -v vcf -p poplist [-w window [default:100000]] [-s step [default:50000]] [-o output [default:Dxy.txt]]
# script 
# module load python/3.6.4
# for i in {1..10}; do
#   python3 Dxy_calculate.py -v 'SL_OBOYS2_chr'$i'.recode.vcf.gz' -p SL_OBOYS2.txt -w 500000 -s 50000 -o 'SL_OBOYS2_chr'$i'.dxy'
#   python3 Dxy_calculate.py -v 'CS_DEBY_chr'$i'.recode.vcf.gz' -p CS_DEBY.txt -w 500000 -s 50000 -o 'CS_DEBY_chr'$i'.dxy'
#   python3 Dxy_calculate.py -v 'CS_NEH_chr'$i'.recode.vcf.gz' -p CS_NEH.txt -w 500000 -s 50000 -o 'CS_NEH_chr'$i'.dxy'
# done
# load the R function
setwd("~/Documents/Ryan_workplace/CVreseq_hudson_fst_nobiginvers/Dxy_plot")
source("manhattan.R")
library(export)

#plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_noinvers.1000bp.s200.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
#mid_pos <- round((DT1$window_start + DT1$window_end)/2)
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_NEH_noinvers.1000bp.s200.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
#mid_pos <- round((DT2$window_start + DT2$window_end)/2)
mid_pos <- DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "CS_UMFS_noinvers.1000bp.s200.csv"
DT3 = read.delim(name3, header = TRUE, sep=',')
#mid_pos <- round((DT3$window_start + DT3$window_end)/2)
mid_pos <- DT3$mid
id3 = paste0(DT3$scaffold,'_',mid_pos)
DT3 <- as.data.frame(cbind(DT3,mid_pos,id3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]

common_id = intersect(intersect(DT1$id1,DT2$id2),DT3$id3)
common_id <- as.vector(common_id)
sub_DT1 = DT1[DT1$id1 %in% common_id,]
sub_DT2 = DT2[DT2$id2 %in% common_id,]
sub_DT3 = DT3[DT3$id3 %in% common_id,]

sub_DT1 = sub_DT1[match(common_id, sub_DT1$id1),]
sub_DT2 = sub_DT2[match(common_id, sub_DT2$id2),]
sub_DT3 = sub_DT3[match(common_id, sub_DT3$id3),]

dat <- data.frame(chr=sub_DT1$scaffold, start=sub_DT1$start, start=sub_DT1$end, mid_pos=sub_DT1$mid, SNP=common_id, D1_dxy=sub_DT1$dxy_CS_DEBY, D2_dxy=sub_DT2$dxy_CS_NEH, D3_dxy=sub_DT3$dxy_CS_UMFS)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D1_dxy <- as.numeric(dat$D1_dxy)
dat$D2_dxy <- as.numeric(dat$D2_dxy)
dat$D3_dxy <- as.numeric(dat$D3_dxy)
write.table(dat, file = "10K.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')
###################### start dxy plot ######################
jpeg("Dxy_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="mid_pos",p="D1_dxy", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.015), xlim = c(32800000, 32950000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-DEBY Dxy (10000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D2_dxy", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.015), xlim = c(32800000, 32950000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-NEH Dxy (10000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_dxy", subset(dat, chr == 5),  highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.015), xlim = c(32800000, 32950000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", cex.lab=1.5, main = "CS_UMFS Dxy (10000 bp/window)", cex.main=1.5)
#graph2ppt(file="Dxy_1K.pptx", width=9, height=6.5)
dev.off()
#  print("Pi Plotting done")
#}



#plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_noinvers.1000bp.s200.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
#mid_pos <- round((DT1$window_start + DT1$window_end)/2)
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
p1 = DT1$pi_CS - DT1$pi_DEBY
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1, p1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_NEH_noinvers.1000bp.s200.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
#mid_pos <- round((DT2$window_start + DT2$window_end)/2)
mid_pos <- DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
p2 = DT2$pi_CS - DT2$pi_NEH
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2,p2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "CS_UMFS_noinvers.1000bp.s200.csv"
DT3 = read.delim(name3, header = TRUE, sep=',')
#mid_pos <- round((DT3$window_start + DT3$window_end)/2)
mid_pos <- DT3$mid
id3 = paste0(DT3$scaffold,'_',mid_pos)
p3 = DT3$pi_CS - DT3$pi_UMFS
DT3 <- as.data.frame(cbind(DT3,mid_pos,id3,p3))
# remove NA rows in the dxy column
DT3 <- DT3[complete.cases(DT3), ]

common_id = intersect(intersect(DT1$id1,DT2$id2),DT3$id3)
common_id <- as.vector(common_id)
sub_DT1 = DT1[DT1$id1 %in% common_id,]
sub_DT2 = DT2[DT2$id2 %in% common_id,]
sub_DT3 = DT3[DT3$id3 %in% common_id,]

sub_DT1 = sub_DT1[match(common_id, sub_DT1$id1),]
sub_DT2 = sub_DT2[match(common_id, sub_DT2$id2),]
sub_DT3 = sub_DT3[match(common_id, sub_DT3$id3),]

dat <- data.frame(chr=sub_DT1$scaffold, start=sub_DT1$start, start=sub_DT1$end, mid_pos=sub_DT1$mid, SNP=common_id, D1_p=sub_DT1$p1, D2_p=sub_DT2$p2, D3_p=sub_DT3$p3)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D1_p <- as.numeric(dat$D1_p)
dat$D2_p <- as.numeric(dat$D2_p)
dat$D3_p <- as.numeric(dat$D3_p)
write.table(dat, file = "1K.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')
###################### start dxy plot ######################
jpeg("pi_1K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="mid_pos",p="D1_p", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(-0.005, 0.005), xlim = c(32900000, 32920000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="pi diff", xlab="", cex.lab=1.5, main = "CS-DEBY pi diff (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D2_p", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(-0.005, 0.005), xlim = c(32900000, 32920000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="pi diff", xlab="", cex.lab=1.5, main = "CS-NEH pi diff (1000 bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_p", subset(dat, chr == 5),  highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(-0.005, 0.005), xlim = c(32900000, 32920000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="pi diff", cex.lab=1.5, main = "CS-UMFS pi diff (1000 bp/window)", cex.main=1.5)
#graph2ppt(file="Dxy_1K.pptx", width=9, height=6.5)
dev.off()
#  print("Pi Plotting done")
#}

