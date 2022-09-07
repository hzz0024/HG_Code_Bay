rm(list=ls())


################################################################
############# find PCAdapt and outlier shared results ##########
################################################################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/09_shared_outliers/pop_477/")
win_s = 5000
outflank <- read.delim(paste0("pop_n_477_outflank_BP_q05_n_276.bed"), header = FALSE, sep='\t')
outflank_df <- data.frame(outflank$V1, outflank$V2-win_s, outflank$V3+win_s)
write.table(outflank_df, file = "pop_n_477_outflank_BP_q05_n_276_1k.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

PCAdapt <- read.delim(paste0("pop_n_477_PCAdapt_BP_q05_n_386.bed"), header = FALSE, sep='\t')
PCAdapt_df <- data.frame(PCAdapt$V1, PCAdapt$V2-win_s, PCAdapt$V3+win_s)
write.table(PCAdapt_df, file = "pop_n_477_PCAdapt_BP_q05_n_386_1k.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# number of shared SNPs
length(intersect(outflank$V4, PCAdapt$V4))

# headname1 = 'pop_n_477_outflank_BP_q05_n_276_1k.bed'
# headname2 = 'pop_n_477_PCAdapt_BP_q05_n_386_1k.bed'
# DT1 = read.delim(headname1, header = FALSE, sep='\t')
# DT2 = read.delim(headname2, header = FALSE, sep='\t')
# colnames(DT1) = c("Chromosome", "Start", "End")
# colnames(DT2) = c("Chromosome", "Start", "End")
# ST = c()
# ED = c()
# CHR = c()
# for(i in seq(length(DT1[,1]))){
#   st = DT1$Start[i]
#   ed = DT1$End[i]
#   chr = DT1$Chromosome[i]
#   for(j in seq(length(DT2[,1]))){
#     st_ = DT2$Start[j]
#     ed_ = DT2$End[j]
#     chr_ = DT2$Chromosome[j]
#     if( (ed>=st_&&st<=ed_) && chr==chr_){
#       ST = c(ST, min(st,st_))
#       ED = c(ED,  max(ed,ed_))
#       CHR = c(CHR, chr)
#     }
#   }
# }
# overlaps = data.frame(chr=CHR, start=ST, end=ED, Feature="shared_outlier", ZFst=10)
# write.table(overlaps, file = paste0(headname1,"_",headname2,".outlier.igv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
# for(i in seq(10)) 
#   overlaps$chr[overlaps$chr==i] = chr_str_list[i]
# overlap_dt <- overlaps[with(overlaps, order(chr, start)),]
# write.table(overlap_dt, file = "pop_n_477_pcadapt_outflank.shared.10k.outlier.igv", quote = FALSE, row.names = FALSE, col.names = TRUE)

# overlap_snp <- intersect(outflank$V4, PCAdapt$V4)
# overlap_dt <- outflank[which(outflank$V4 %in% overlap_snp),]

shared_SNPs <- intersect(outflank$V4, PCAdapt$V4)
outflank_PCAdapt_shared_SNPs = PCAdapt[which(PCAdapt$V4 %in% shared_SNPs),]

chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(10)) 
  outflank_PCAdapt_shared_SNPs$V1[outflank_PCAdapt_shared_SNPs$V1==i] = chr_str_list[i]
overlap_dt <- outflank_PCAdapt_shared_SNPs[with(outflank_PCAdapt_shared_SNPs, order(V1, V2)),]
write.table(overlap_dt, file = "pop_n_477_pcadapt_outflank.shared.outlier.igv", quote = FALSE, row.names = FALSE, col.names = TRUE)

############# find shared outliers from merged results ##########
#################################################################
library(hash)
options(scipen=999)
library(gtools)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/09_shared_outliers/pop_477")
headname1 = 'pop_n_477_pcadapt_outflank.shared.outlier.igv'
headname2 = 'Dom_Wild.sliding.zfst.outlier.merged.igv'
DT1 = read.delim(headname1, header = TRUE, sep=' ')
colnames(DT1) = c("Chromosome", "Start", "End", "ID")
DT2 = read.delim(headname2, header = TRUE, sep='\t')
ST = c()
ED = c()
CHR = c()
for(i in seq(length(DT1[,1]))){
  st = DT1$Start[i]
  ed = DT1$End[i]
  chr = DT1$Chromosome[i]
  for(j in seq(length(DT2[,1]))){
    st_ = DT2$Start[j]
    ed_ = DT2$End[j]
    chr_ = DT2$Chromosome[j]
    if( (ed>=st_&&st<=ed_) && chr==chr_){
      ST = c(ST, min(st,st_))
      ED = c(ED,  max(ed,ed_))
      CHR = c(CHR, chr)
    }
  }
}
#print(CHR)
#print(ST)
#print(ED)
if(length(CHR) > 0){
  overlaps = data.frame(chr=CHR, start=ST, end=ED, Feature="Merged_outlier", ZFst=10)
  overlaps = overlaps[order(overlaps$chr, overlaps$start),]
  colnames(overlaps)=c('Chromosome', 'Start', 'End', 'Feature', paste0(headname1,"-",headname2,'_ZFst'))
  overlaps <- unique(overlaps)
  print(length(overlaps$Chromosome))
  write.table(overlaps, file = paste0(headname1,"_",headname2,".outlier.igv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

bedtools = "/Users/HG/Downloads/bedtools2/bin/bedtools";
system(paste(bedtools, " merge -i ",paste0(headname1,"_",headname2,".outlier.igv"), " -c 1 -o count > ",paste0(headname1,"_",headname2,".outlier.merged.igv"),  sep=""))
#bedtools merge -i pop_n_282_shared.outlier.merged.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.igv -c 1 -o count > pop_n_282_shared.outlier.merged.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv




