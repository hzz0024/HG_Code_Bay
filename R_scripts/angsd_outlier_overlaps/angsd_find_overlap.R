###############################################################
### Bedtools : find the intersection between SNPs and genes ###
###############################################################

# the code is adopted from the tutorial here https://github.com/clairemerot/physalia_adaptation_course/tree/2021/05_day5

###########################################
### format outlier files into bed files ###
###########################################


#######################
#  make bed file      #
# bed format files are used to find the intrersection regions from different contrasts
#load file
setwd("~/Documents/Ryan_workplace/DelBay_adult/11_SGS")

bed_format <- function(input, output){
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

bed_format("Del19_FDR_outlier.list", "./bed/Del19_FDR_10K.bed" )
bed_format("Del20_FDR_outlier.list", "./bed/Del20_FDR_10K.bed" )
bed_format("NB_HC_FDR_outlier.list", "./bed/NB_HC_FDR_10K.bed" )
bed_format("SR_HC_FDR_outlier.list", "./bed/SR_HC_FDR_10K.bed" )

# do the same thing for whole SNPs

outlier_temp_SGS<-read.table("ps_Del19_challenge.txt", header=F)
#what's its dimension?
dim(outlier_temp_SGS)
head(outlier_temp_SGS)
#which size around the SNP
window<-10000
#add a vector with start position
outlier_temp_SGS$start<-outlier_temp_SGS$V2-(window/2)
#oups it can't be ngative! replace negative by 0
outlier_temp_SGS$start[outlier_temp_SGS$start<0]<-0 
#add a vector with stop position
outlier_temp_SGS$stop<-outlier_temp_SGS$V2+(window/2)
outlier_temp_SGS$id_snp <- paste0(outlier_temp_SGS$V1,'_',outlier_temp_SGS$V2)
#have a look
head(outlier_temp_SGS)
#which columns shoud we keep?
outlier_temp_SGS_bed<-outlier_temp_SGS[,c(1,8,9,10)]
colnames(outlier_temp_SGS_bed) <- c("chromosome", "start", "stop", "id_snp")
head(outlier_temp_SGS_bed)
#save your file
write.table(outlier_temp_SGS_bed, "./bed/all_snps.bed", row.names=F, sep="\t", quote=F,col.names=F)

############# find shared outliers from non-merged results ############# 
library(hash)
options(scipen=999)
library(gtools)

#heads = c("SL_OBOYS2","SL_LOLA","NEH_UMFS","CS_UMFS",'CS_DEBY', 'CS_NEH','CL_OBOYS2',"CS_HC")
#heads = c("CS_HC",'HC_CLP', 'HCVA_CLP','CS_HCVA')
heads = c("CS_DEBY",'CS_NEH', 'CS_UMFS','CL_OBOYS2')
combs = combinations(4, 2, 1:4)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0(headname1,"_noinvers.sliding.zfst.outlier.igv"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0(headname2,"_noinvers.sliding.zfst.outlier.igv"), header = TRUE, sep='\t')
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
      if( (st == st_ && ed == ed_) && chr==chr_){
        ST = c(ST, st)
        ED = c(ED, ed)
        CHR = c(CHR, chr)
      }
    }
  }
  if(length(CHR) > 0){
    overlaps = data.frame(chr=CHR, start=ST, end=ED, Feature="Shared_outlier", ZFst=10)
    overlaps = overlaps[order(overlaps$chr, overlaps$start),]
    colnames(overlaps)=c('Chromosome', 'Start', 'End', 'Feature', paste0(headname1,"-",headname2,'_ZFst'))
    print(length(overlaps$Chromosome))
    write.table(overlaps, file = paste0(headname1,"-",headname2,".outlier.shared.igv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}
############# find shared outliers from merged results ############# 
library(hash)
options(scipen=999)
library(gtools)
setwd("~/Documents/Ryan_workplace/CVreseq_hudson_fst_nobiginvers/Zfst_1K")
#heads = c("SL_OBOYS2","SL_LOLA","NEH_UMFS","CS_UMFS",'CS_DEBY', 'CS_NEH','CL_OBOYS2',"CS_HC")
#heads = c("CS_DEBY",'CS_NEH', 'CS_UMFS','CL_OBOYS2')
heads = c("CS_HC", 'HCVA_CLP', "HI_SM")
combs = combinations(3, 2, 1:3)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0(headname1,"_noinvers.sliding.zfst.outlier.merged.igv"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0(headname2,"_noinvers.sliding.zfst.outlier.merged.igv"), header = TRUE, sep='\t')
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
        ST = c(ST, max(st,st_))
        ED = c(ED,  min(ed,ed_))
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
    print(length(overlaps$Chromosome))
    write.table(overlaps, file = paste0(headname1,"-",headname2,".outlier.merged.igv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

name3 = "CS_UMFS_noinvers.sliding.zfst.outlier.igv"
DT3 = read.delim(name3, header = TRUE, sep='\t')

ST = c()
ED = c()
CHR = c()
for(i in seq(length(DT3[,1]))){
  st = DT3$Start[i]
  ed = DT3$End[i]
  chr = DT3$Chromosome[i]
  for(j in seq(length(overlaps[,1]))){
    st_ = overlaps$start[j]
    ed_ = overlaps$end[j]
    chr_ = overlaps$chr[j]
    if( (ed>=st_&&st<=ed_) && chr==chr_){
      ST = c(ST, max(st,st_))
      ED = c(ED,  min(ed,ed_))
      CHR = c(CHR, chr)
    }
  }
}
overlaps = data.frame(chr=CHR, start=ST, end=ED)
overlaps = overlaps[order(overlaps$chr, overlaps$start),]

################## find shared ouliter and merging (wrong) ################## 
#setwd("~/Documents/Ryan_workplace/CVreseq_hudson_fst_nobiginvers/Zfst_1K")
# 
# shared_outlier <- function(headname1,headname2){
#   #headname1 = "CS_DEBY"
#   DT1 = read.delim(paste0(headname1,"_noinvers.sliding.zfst.outlier.igv"), header = TRUE, sep='\t')
#   DT1$id1 = paste0(DT1$Chromosome,'_',DT1$Start,'_',DT1$End)
#   #headname2 = "CS_NEH"
#   DT2 = read.delim(paste0(headname2,"_noinvers.sliding.zfst.outlier.igv"), header = TRUE, sep='\t')
#   DT2$id2 = paste0(DT2$Chromosome,'_',DT2$Start,'_',DT2$End)
#   common = intersect(DT1$id1, DT2$id2)
#   #intersect(DT1$id1, intersect(DT2$id2, DT3$id3))
#   length(common)
#   shared = DT1[DT1$id1 %in% common,]
#   ZFst2 = DT2[DT2$id2 %in% common,][5]
#   share_DT = data.frame(Chromosome=shared$Chromosome, Start=shared$Start, End=shared$End, Feature="shared", ZFst1=shared[5], ZFst2=ZFst2)
#   colnames(share_DT)=c('Chromosome', 'Start', 'End', 'Feature', 'ZFst1','ZFst2')
#   #====================MERGE======================
#   #outlier_DT = share_DT[,-c(4,5,6)] #delete two columns
#   outlier_DT = share_DT[order(share_DT$Chromosome, share_DT$Start),]
#   library(hash)
#   dict = hash()
#   Zdict1 = hash()
#   Zdict2 = hash()
#   #for(chr in unique(outlier_DT$Chromosome)){
#   #  Zdict[[chr]] = c()
#   #}
#   zFst_set1 = c()
#   zFst_set2 = c()
#   for(i in seq(length(outlier_DT[,1]))){
#     chr = outlier_DT$Chromosome[i]
#     st = outlier_DT$Start[i]
#     ed = outlier_DT$End[i]
#     ranges = dict[[chr]]
#     if(length(ranges)>0 && ranges[[length(ranges)]][2]>=st){ #overlap
#       last_range = ranges[[length(ranges)]]
#       ranges = ranges[-length(ranges)] #remove last one
#       dict[[chr]] = append(ranges, list(c(last_range[1], max(last_range[2], ed))))
#       # previous zFst set add one
#       zFst_set1 = c(zFst_set1, outlier_DT$ZFst1[i])
#       zFst_set2 = c(zFst_set2, outlier_DT$ZFst2[i])
#       #print(zFst_set)
#     }else{
#       dict[[chr]] = append(ranges, list(c(st, ed)))
#       #previous zFst set calculate median, write,. reset zFst set
#       if(i==1){
#         zFst_set1 = c(zFst_set1, outlier_DT$ZFst1[i])
#         zFst_set2 = c(zFst_set2, outlier_DT$ZFst2[i])
#       }else{
#         max_zFst1 = max(zFst_set1)
#         max_zFst2 = max(zFst_set2)
#         Zdict1[[outlier_DT$Chromosome[i-1]]] = c(Zdict1[[outlier_DT$Chromosome[i-1]]], max_zFst1)
#         Zdict2[[outlier_DT$Chromosome[i-1]]] = c(Zdict2[[outlier_DT$Chromosome[i-1]]], max_zFst2)
#         zFst_set1 = c(outlier_DT$ZFst1[i])
#         zFst_set2 = c(outlier_DT$ZFst2[i])
#       }
#     }
#   }
#   max_zFst1 = max(zFst_set1)
#   max_zFst2 = max(zFst_set2)
#   Zdict1[[chr]] = c(Zdict1[[chr]], max_zFst1)
#   Zdict2[[chr]] = c(Zdict2[[chr]], max_zFst2)
#   #===================================================
#   CHR = c()
#   ST = c()
#   ED = c()
#   ZF1 = c()
#   ZF2 = c()
#   for(chr in unique(outlier_DT$Chromosome)){
#     ranges = dict[[chr]] 
#     zFsts1 = Zdict1[[chr]]
#     zFsts2 = Zdict2[[chr]]
#     for(i in seq(length(ranges))){
#       CHR = c(CHR, chr)
#       ST = c(ST, ranges[[i]][1])
#       ED = c(ED, ranges[[i]][2])
#       ZF1 = c(ZF1, zFsts1[[i]])
#       ZF2 = c(ZF2, zFsts2[[i]])
#     }
#   }
#   outlier_DT = data.frame(Chromosome=CHR, Start=ST, End=ED, Feature="merged",ZFst1=ZF1,ZFst2=ZF2)
#   print(paste0("Number of merged outlier is ", length(outlier_DT$Feature)))
#   write.table(outlier_DT, file = paste0(headname1, "_", headname2, ".zfst.shared.merged.igv"), sep = "\t", quote = FALSE,
#               row.names = FALSE, col.names = TRUE)
# }
# 
# shared_outlier("CS_DEBY", "CS_NEH")
# shared_outlier("CS_DEBY", "CS_UMFS")




