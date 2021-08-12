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
name1 = "CS_DEBY_noinvers.10000bp.s2000.csv"
DT1 = read.delim(name1, header = TRUE, sep=',')
#mid_pos <- round((DT1$window_start + DT1$window_end)/2)
mid_pos <- DT1$mid
id1 = paste0(DT1$scaffold,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_NEH_noinvers.10000bp.s2000.csv"
DT2 = read.delim(name2, header = TRUE, sep=',')
#mid_pos <- round((DT2$window_start + DT2$window_end)/2)
mid_pos <- DT2$mid
id2 = paste0(DT2$scaffold,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "CS_UMFS_noinvers.10000bp.s2000.csv"
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
########### chr 3 #############
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
########### chr 4 #############
#plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_chr4.dxy"
DT1 = read.delim(name1, header = TRUE, sep='\t')
mid_pos <- round((DT1$window_start + DT1$window_end)/2)
id1 = paste0(DT1$chrom,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_NEH_chr4.dxy"
DT2 = read.delim(name2, header = TRUE, sep='\t')
mid_pos <- round((DT2$window_start + DT2$window_end)/2)
id2 = paste0(DT2$chrom,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "SL_OBOYS2_chr4.dxy"
DT3 = read.delim(name3, header = TRUE, sep='\t')
mid_pos <- round((DT3$window_start + DT3$window_end)/2)
id3 = paste0(DT3$chrom,'_',mid_pos)
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

dat <- data.frame(chr=sub_DT1$chrom, start=sub_DT1$window_start, start=sub_DT1$window_end, mid_pos=sub_DT1$mid_pos, SNP=common_id, D1_dxy=sub_DT1$Dxy, D2_dxy=sub_DT2$Dxy, D3_dxy=sub_DT3$Dxy)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D1_dxy <- as.numeric(dat$D1_dxy)
dat$D2_dxy <- as.numeric(dat$D2_dxy)
dat$D3_dxy <- as.numeric(dat$D3_dxy)
write.table(dat, file = "500K_chr4.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest_chr4.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')
###################### start dxy plot ######################
#jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="mid_pos",p="D1_dxy", subset(dat, chr == 4), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-DEBY Dxy (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D2_dxy", subset(dat, chr == 4), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-NEH Dxy (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_dxy", subset(dat, chr == 4),  highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", cex.lab=1.5, main = "SL-OBOYS Dxy (500K bp/window)", cex.main=1.5)
graph2ppt(file="Dxy_chr4_500K.pptx", width=9, height=6.5)
#dev.off()
#  print("Pi Plotting done")
#}

########### chr 5 #############
#plot_pi <- function(name1, name2, plot_name, pop){
# load the original windowed pi
name1 = "CS_DEBY_chr5.dxy"
DT1 = read.delim(name1, header = TRUE, sep='\t')
mid_pos <- round((DT1$window_start + DT1$window_end)/2)
id1 = paste0(DT1$chrom,'_',mid_pos)
DT1 <- as.data.frame(cbind(DT1,mid_pos, id1))
# remove NA rows in the dxy column
DT1 <- DT1[complete.cases(DT1), ]

name2 = "CS_NEH_chr5.dxy"
DT2 = read.delim(name2, header = TRUE, sep='\t')
mid_pos <- round((DT2$window_start + DT2$window_end)/2)
id2 = paste0(DT2$chrom,'_',mid_pos)
DT2 <- as.data.frame(cbind(DT2,mid_pos,id2))
# remove NA rows in the dxy column
DT2 <- DT2[complete.cases(DT2), ]

name3 = "SL_OBOYS2_chr5.dxy"
DT3 = read.delim(name3, header = TRUE, sep='\t')
mid_pos <- round((DT3$window_start + DT3$window_end)/2)
id3 = paste0(DT3$chrom,'_',mid_pos)
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

dat <- data.frame(chr=sub_DT1$chrom, start=sub_DT1$window_start, start=sub_DT1$window_end, mid_pos=sub_DT1$mid_pos, SNP=common_id, D1_dxy=sub_DT1$Dxy, D2_dxy=sub_DT2$Dxy, D3_dxy=sub_DT3$Dxy)
dat$chr <- as.numeric(dat$chr)
dat$mid_pos <- as.numeric(dat$mid_pos)
dat$D1_dxy <- as.numeric(dat$D1_dxy)
dat$D2_dxy <- as.numeric(dat$D2_dxy)
dat$D3_dxy <- as.numeric(dat$D3_dxy)
write.table(dat, file = "500K_chr5.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
# import the highlight snps list
file1 = 'snpsOfInterest_chr5.csv' 
h1 = read.delim(file1, header = FALSE, sep=',')
h1 = as.list(h1)
print('Prepare Data Done')
###################### start dxy plot ######################
#jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1))
manhattan(chr="chr",bp="mid_pos",p="D1_dxy", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-DEBY Dxy (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D2_dxy", subset(dat, chr == 5), highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", xlab="", cex.lab=1.5, main = "CS-NEH Dxy (500K bp/window)", cex.main=1.5)
manhattan(chr="chr",bp="mid_pos",p="D3_dxy", subset(dat, chr == 5),  highlight1 = h1$V1, logp=FALSE, cex.axis = 1, ylim = c(0, 0.0005), #xlim = c(50000000, 70000000),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="Dxy", cex.lab=1.5, main = "SL-OBOYS Dxy (500K bp/window)", cex.main=1.5)
graph2ppt(file="Dxy_chr5_500K.pptx", width=9, height=6.5)
#dev.off()
#  print("Pi Plotting done")
#}


