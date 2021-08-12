
setwd("~/Documents/Ryan_workplace/CVreseq_hudson_fst_nobiginvers")
library(hash)
source("manhattan.R")
library(hash)
options(scipen=999)

############################## independent function for sliding window ##############################
############################## independent function for sliding window ##############################
############################## independent function for sliding window ##############################
# headname = "SL_OBOYS2_noinvers."
# titlename = 'SL vs OBOYS2'
hudson_fst <- function(headname, titlename){
  jpeg(paste0(headname,"sliding.zfst.hudson.jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(1,1))
  for(win in c(1000)){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    DT[,9][DT[,9]<0] = 0.000001 #@chnage
    zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE))# @change
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9], zfst = zfst) # @change
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    dat$zfst <- as.numeric(dat$zfst)
    # calculate the 99.9% quantile
    thred=quantile(dat$zfst, 0.999)
    cnt = length(dat[dat$zfst>thred[[1]],1])
    outlier = dat[dat$zfst>thred[[1]],5]
    #cnt = sum(dat$zfst>5)
    #outlier = dat[dat$zfst>5,5]
    print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))
    manhattan(dat, chr="chr",bp="mid_pos",p="zfst", highlight1 = outlier, logp=FALSE, cex.axis = 1, ylim = c(0, max(dat$zfst)+0.2), #subset(dat, chr == 8)
              col=c("grey","black"),genomewideline=F, suggestiveline=F,
              ylab="ZFst", cex.lab=1.5, main = paste0(titlename, " ZFst window size:", win, ' bp'), cex.main=1.5)
    dev.off()
  }
  
  chrom_hash <- hash()
  for(chr in seq(10)){
    key = as.character(chr)
    chrom_hash[[key]] = list()
  }
  windows = c(1000)
  for(win in windows){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    DT[,9][DT[,9]<0] = 0.000001
    zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE)) #@change
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9], zfst = zfst)
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    dat$zfst <- as.numeric(dat$zfst)
    thred=quantile(dat$zfst, 0.999)
    cnt = length(dat[dat$zfst>thred[[1]],1])
    outlier = dat[dat$zfst>thred[[1]],5]
    zfsts = dat[dat$zfst>thred[[1]],7]
    #cnt = sum(dat$zfst>5)
    #outlier = dat[dat$zfst>5,5]
    #zfsts = dat[dat$zfst>5,7]
    #print(outlier)
    ss <- strsplit(outlier, '_')
    for(i in seq(length(ss))){
      s = ss[[i]]
      chr_ = s[1]
      mid_ = as.integer(s[2])
      range = c(mid_-win/2, mid_+win/2, zfsts[i])
      i = i + 1
      pre_ranges = chrom_hash[[chr_]]
      if(win==windows[1]){
        # init
        chrom_hash[[chr_]] = append(chrom_hash[[chr_]], list(range))
      }else{
        OVERLAP = FALSE
        for(small_range in pre_ranges){
          if( range[2]>=small_range[1] && range[1]<=small_range[2] ){
            #if( (range[1]>=small_range[1] && range[1]<=small_range[2]) || ( range[2]>=small_range[1] && range[2]<=small_range[2]) ){
            OVERLAP = TRUE
            break
          }
        }
        if(!OVERLAP){
          chrom_hash[[chr_]] = append(chrom_hash[[chr_]], list(range))
        }
      }
    }
    if(win==windows[length(windows)]){
      #print(chrom_hash)
    }
    #print(unlist(new_chrom_hash[['8']]))
  }
  
  # print final ranges to bed file
  fileConn <- file(paste0(headname,'sliding.zfst.outlier.igv'))
  ls = c()
  chrs = c()
  sts = c()
  eds = c()
  ZFSTs = c()
  cnt = 0
  for(chr in seq(10)){
    chr = as.character(chr)
    if(length(chrom_hash[[chr]])>0){
      ranges = chrom_hash[[chr]]
      for(range in ranges){
        cnt = cnt + 1
        l = paste0(chr, '\t', range[1], '\t', range[2])
        chrs = c(chrs, as.integer(chr))
        sts = c(sts, range[1])
        eds = c(eds, range[2])
        ZFSTs = c(ZFSTs, range[3])
        ls = c(ls, l)
      }        
    }
  }
  #print(cnt)
  writeLines(ls, fileConn)
  close(fileConn)
  
  outlier_DT = data.frame(chr=chrs, st=sts, ed=eds, SNP=paste0(chrs,"_",round((sts + eds)/2)), zfst=ZFSTs)
  outlier_DT = outlier_DT[order(outlier_DT$chr, outlier_DT$st),]
  print(paste0("Number of outliers in ", titlename, " is ", length(outlier_DT$SNP)))
  #jpeg("CS_NEH_noinvers_overlaps.jpg", width = 16, height = 9, units = 'in', res = 300)
  win = 1000
  name = paste0(headname, win, "bp.", "s", win/5, ".csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  DT[,9][DT[,9]<0] = 0.000001 #@chnage
  DT$zfst <- (DT[,9] - mean(DT[,9]))/(sd(DT[,9], na.rm = FALSE))
  dat <- data.frame(chr=c(DT$scaffold, outlier_DT$chr), start=c(DT$start, outlier_DT$st), 
                    end=c(DT$end, outlier_DT$ed), mid_pos=c(DT$mid_pos, round((outlier_DT$st + outlier_DT$ed)/2)),  
                    zfst=c(DT[,12], outlier_DT$zfst), SNP=c(DT$id, paste0(outlier_DT$SNP,'_outlier'))) 
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$zfst <- as.numeric(dat$zfst)
  
  #manhattan(dat, chr="chr",bp="mid_pos",p="zfst", highlight1 = paste0(outlier_DT[outlier_DT$ed-outlier_DT$st==1000,4],'_outlier'), highlight2=paste0(outlier_DT[outlier_DT$ed-outlier_DT$st!=1000,4],'_outlier'), logp=FALSE, cex.axis = 1, ylim = c(0, max(dat$zfst)+0.2), #subset(dat, chr == 6), xlim = c(21400000, 21420000),
            #col=c("grey","black"),genomewideline=F, suggestiveline=F,
            #ylab="ZFst", cex.lab=1.5, main = paste0(titlename," union outliers based 1000, 5000, 10000 bp window tests"), cex.main=1.5)
  #dev.off()
  # output bed file
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature', paste0(titlename,'_ZFst'))
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_DT$Chromosome[outlier_DT$Chromosome==i] = chr_str_list[i]
  for(win in windows)
    outlier_DT$Feature[outlier_DT$End - outlier_DT$Start == win] = win
  write.table(outlier_DT, file = paste0(headname,"sliding.zfst.outlier.igv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  #====================MERGE======================
  #outlier_DT = outlier_DT[,-c(4,5)] #delete two columns
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature', 'ZFst')
  outlier_DT = outlier_DT[order(outlier_DT$Chromosome, outlier_DT$Start),]
  library(hash)
  dict = hash()
  Zdict = hash()
  #for(chr in unique(outlier_DT$Chromosome)){
  #  Zdict[[chr]] = c()
  #}
  zFst_set = c()
  for(i in seq(length(outlier_DT[,1]))){
    chr = outlier_DT$Chromosome[i]
    st = outlier_DT$Start[i]
    ed = outlier_DT$End[i]
    ranges = dict[[chr]]
    if(length(ranges)>0 && ranges[[length(ranges)]][2]>=st){ #overlap
      last_range = ranges[[length(ranges)]]
      ranges = ranges[-length(ranges)] #remove last one
      dict[[chr]] = append(ranges, list(c(last_range[1], max(last_range[2], ed))))
      # previous zFst set add one
      zFst_set = c(zFst_set, outlier_DT$ZFst[i])
      #print(zFst_set)
    }else{
      dict[[chr]] = append(ranges, list(c(st, ed)))
      #previous zFst set calculate median, write,. reset zFst set
      if(i==1){
        zFst_set = c(zFst_set, outlier_DT$ZFst[i])
      }else{
        max_zFst = max(zFst_set)
        Zdict[[outlier_DT$Chromosome[i-1]]] = c(Zdict[[outlier_DT$Chromosome[i-1]]], max_zFst)
        zFst_set = c(outlier_DT$ZFst[i])
      }
    }
  }
  max_zFst = max(zFst_set)
  Zdict[[chr]] = c(Zdict[[chr]], max_zFst)
  
  #for(chr in unique(outlier_DT$Chromosome)){
  #print('---------')
  #print(chr)
  #print(Zdict[[chr]])
  #print(length(dict[[chr]]))
  #}
  #===================================================
  CHR = c()
  ST = c()
  ED = c()
  ZF = c()
  for(chr in unique(outlier_DT$Chromosome)){
    ranges = dict[[chr]] 
    zFsts = Zdict[[chr]]
    for(i in seq(length(ranges))){
      CHR = c(CHR, chr)
      ST = c(ST, ranges[[i]][1])
      ED = c(ED, ranges[[i]][2])
      ZF = c(ZF, zFsts[[i]])
    }
  }
  outlier_DT = data.frame(Chromosome=CHR, Start=ST, End=ED, Feature="merged",ZFst=ZF)
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature', paste0(titlename,'_ZFst'))
  write.table(outlier_DT, file = paste0(headname,"sliding.zfst.outlier.merged.igv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}


############################## formal run ##############################
############################## formal run ##############################
############################## formal run ##############################

hudson_fst("SL_OBOYS2_noinvers.", "SL_OBOYS2")
hudson_fst("SL_LOLA_noinvers.", "SL_LOLA")
hudson_fst("NEH_UMFS_noinvers.", "NEH_UMFS")
hudson_fst("CS_UMFS_noinvers.", "CS_UMFS")
hudson_fst("CS_NEH_noinvers.", "CS_NEH")
hudson_fst("CS_DEBY_noinvers.", "CS_DEBY")
hudson_fst("CL_OBOYS2_noinvers.", "CL_OBOYS2")
hudson_fst("CS_SL_noinvers.", "CS_SL")

hudson_fst("CS_HC_noinvers.", "CS_HC")
hudson_fst("HC_CLP_noinvers.", "HC_CLP")
hudson_fst("HCVA_CLP_noinvers.", "HCVA_CLP")
hudson_fst("CS_HCVA_noinvers.", "CS_HCVA")
hudson_fst("HI_SM_noinvers.", "HI_SM")

