setwd("/Volumes/cornell/CVreseq_hudson_fst")
library(hash)
source("manhattan.R")
library(hash)
install.packages("KRIS")
options(scipen=999)
#mean1 <- function(x) (x[1]+x[2])/2
#all_outlier = c()
#for(chr in seq(10)){
#  chr = as.character(chr)
#  if(length(chrom_hash[[chr]])>0){
#    overlap_outlier = unlist(lapply(chrom_hash[[chr]], mean1))
#    overlap_outlier = paste0(chr, '_', overlap_outlier)
#    all_outlier = c(all_outlier, overlap_outlier)
#  }
#}

############################## independent function ##############################
############################## independent function ##############################
############################## independent function ##############################

hudson_fst <- function(headname, titlename){
  jpeg(paste0(headname,"fst.hudson.jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(2,2))
  for(win in c(1000, 5000, 10000)){
    name = paste0(headname, win, "bp.csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9])
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    # calculate the 99.9% quantile
    thred=quantile(dat$fst, 0.999)
    cnt = length(dat[dat$fst>thred[[1]],1])
    outlier = dat[dat$fst>thred[[1]],5]
    print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))
    manhattan(dat, chr="chr",bp="mid_pos",p="fst", highlight1 = outlier, logp=FALSE, cex.axis = 1, ylim = c(0, 1.02), #subset(dat, chr == 8)
              col=c("grey","black"),genomewideline=F, suggestiveline=F,
              ylab="Fst", cex.lab=1.5, main = paste0(titlename, " Fst window size:", win, ' bp'), cex.main=1.5)
  }
  
  chrom_hash <- hash()
  for(chr in seq(10)){
    key = as.character(chr)
    chrom_hash[[key]] = list()
  }
  windows = c(1000, 5000, 10000)
  for(win in windows){
    name = paste0(headname, win, "bp.csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9]) 
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    thred=quantile(dat$fst, 0.999)
    cnt = length(dat[dat$fst>thred[[1]],1])
    outlier = dat[dat$fst>thred[[1]],5]
    fsts = dat[dat$fst>thred[[1]],6]
    #print(outlier)
    ss <- strsplit(outlier, '_')
    for(i in seq(length(ss))){
      s = ss[[i]]
      chr_ = s[1]
      mid_ = as.integer(s[2])
      range = c(mid_-win/2, mid_+win/2, fsts[i])
      i = i + 1
      pre_ranges = chrom_hash[[chr_]]
      if(win==windows[1]){
        # init
        chrom_hash[[chr_]] = append(chrom_hash[[chr_]], list(range))
      }else{
        OVERLAP = FALSE
        for(small_range in pre_ranges){
          #if( range[1]>=big_range[1] && range[2]<=big_range[2] ){
          if( (range[1]>=small_range[1] && range[1]<=small_range[2]) || ( range[2]>=small_range[1] && range[2]<=small_range[2]) ){
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
  
  #print final ranges to bed file
  fileConn <- file(paste0(headname,'final.outlier.igv'))
  ls = c()
  chrs = c()
  sts = c()
  eds = c()
  FSTs = c()
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
        FSTs = c(FSTs, range[3])
        ls = c(ls, l)
      }
    }
  }
  #print(cnt)
  writeLines(ls, fileConn)
  close(fileConn)
  
  outlier_DT = data.frame(chr=chrs, st=sts, ed=eds, SNP=paste0(chrs,"_",round((sts + eds)/2)), fst=FSTs)
  outlier_DT = outlier_DT[order(outlier_DT$chr, outlier_DT$st),]
  
  print(paste0("Number of outliers in ", titlename, " is ", length(outlier_DT$SNP)))
  win = 1000
  name = paste0(headname, win, "bp.csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  
  dat <- data.frame(chr=c(DT$scaffold, outlier_DT$chr), start=c(DT$start, outlier_DT$st), end=c(DT$end, outlier_DT$ed), mid_pos=c(DT$mid_pos, round((outlier_DT$st + outlier_DT$ed)/2)),  fst=c(DT[,9], outlier_DT$fst), SNP=c(DT$id, paste0(outlier_DT$SNP,'_outlier'))) 
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$fst <- as.numeric(dat$fst)
  
  outlier_DT0 = outlier_DT
  
  manhattan(dat, chr="chr",bp="mid_pos",p="fst", highlight1 = paste0(outlier_DT[outlier_DT$ed-outlier_DT$st==1000,4],'_outlier'), highlight2=paste0(outlier_DT[outlier_DT$ed-outlier_DT$st!=1000,4],'_outlier'), logp=FALSE, cex.axis = 1, ylim = c(0, 1.02), #subset(dat, chr == 6), xlim = c(21400000, 21420000),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0(titlename," union outliers based 1000, 5000, 10000 bp window tests"), cex.main=1.5)
  dev.off()
  # output bed file
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature', 'Fst')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_DT$Chromosome[outlier_DT$Chromosome==i] = chr_str_list[i]
  for(win in windows)
    outlier_DT$Feature[outlier_DT$End - outlier_DT$Start == win] = win
  write.table(outlier_DT, file = paste0(headname,"final.outlier.igv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}


############################## formal run ##############################
############################## formal run ##############################
############################## formal run ##############################

hudson_fst("SL_OBOYS2_noinvers.", "SL vs OBOYS2")
hudson_fst("SL_LOLA_noinvers.", "SL vs LOLA")
hudson_fst("NEH_UMFS_noinvers.", "NEH vs UMFS")
hudson_fst("CS_UMFS_noinvers.", "CS vs UMFS")
hudson_fst("CS_NEH_noinvers.", "CS vs NEH")
hudson_fst("CS_HC_noinvers.", "CS vs HC")
hudson_fst("CS_DEBY_noinvers.", "CS vs DEBY")
hudson_fst("CL_OBOYS2_noinvers.", "CL vs OBOYS2")


############################## independent function for sliding window ##############################
############################## independent function for sliding window ##############################
############################## independent function for sliding window ##############################

hudson_fst_sliding <- function(headname, titlename){
  jpeg(paste0(headname,"sliding.fst.hudson.jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(2,2))
  for(win in c(1000, 5000, 10000)){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9])
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    # calculate the 99.9% quantile
    thred=quantile(dat$fst, 0.999)
    cnt = length(dat[dat$fst>thred[[1]],1])
    outlier = dat[dat$fst>thred[[1]],5]
    print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))
    manhattan(dat, chr="chr",bp="mid_pos",p="fst", highlight1 = outlier, logp=FALSE, cex.axis = 1, ylim = c(0, 1.02), #subset(dat, chr == 8)
              col=c("grey","black"),genomewideline=F, suggestiveline=F,
              ylab="Fst", cex.lab=1.5, main = paste0(titlename, " Fst window size:", win, ' bp'), cex.main=1.5)
  }
  
  chrom_hash <- hash()
  for(chr in seq(10)){
    key = as.character(chr)
    chrom_hash[[key]] = list()
  }
  windows = c(1000, 5000, 10000)
  for(win in windows){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    DT <- DT[complete.cases(DT), ]
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9]) 
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    thred=quantile(dat$fst, 0.999)
    cnt = length(dat[dat$fst>thred[[1]],1])
    outlier = dat[dat$fst>thred[[1]],5]
    fsts = dat[dat$fst>thred[[1]],6]
    #print(outlier)
    ss <- strsplit(outlier, '_')
    for(i in seq(length(ss))){
      s = ss[[i]]
      chr_ = s[1]
      mid_ = as.integer(s[2])
      range = c(mid_-win/2, mid_+win/2, fsts[i])
      i = i + 1
      pre_ranges = chrom_hash[[chr_]]
      if(win==windows[1]){
        # init
        chrom_hash[[chr_]] = append(chrom_hash[[chr_]], list(range))
      }else{
        OVERLAP = FALSE
        for(small_range in pre_ranges){
          #if( range[1]>=big_range[1] && range[2]<=big_range[2] ){
          if( (range[1]>=small_range[1] && range[1]<=small_range[2]) || ( range[2]>=small_range[1] && range[2]<=small_range[2]) ){
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
  
  #print final ranges to bed file
  fileConn <- file(paste0(headname,'sliding.final.outlier.igv'))
  ls = c()
  chrs = c()
  sts = c()
  eds = c()
  FSTs = c()
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
        FSTs = c(FSTs, range[3])
        ls = c(ls, l)
      }
    }
  }
  #print(cnt)
  writeLines(ls, fileConn)
  close(fileConn)
  
  outlier_DT = data.frame(chr=chrs, st=sts, ed=eds, SNP=paste0(chrs,"_",round((sts + eds)/2)), fst=FSTs)
  outlier_DT = outlier_DT[order(outlier_DT$chr, outlier_DT$st),]
  
  print(paste0("Number of outliers in ", titlename, " is ", length(outlier_DT$SNP)))
  win = 1000
  name = paste0(headname, win, "bp.", "s", win/5, ".csv")
  DT = read.delim(name, header = TRUE, sep=',')
  mid_pos <- round((DT$start + DT$end)/2)
  id = paste0(DT$scaffold,'_',mid_pos)
  DT <- as.data.frame(cbind(DT,mid_pos, id))
  DT <- DT[complete.cases(DT), ]
  
  dat <- data.frame(chr=c(DT$scaffold, outlier_DT$chr), start=c(DT$start, outlier_DT$st), end=c(DT$end, outlier_DT$ed), mid_pos=c(DT$mid_pos, round((outlier_DT$st + outlier_DT$ed)/2)),  fst=c(DT[,9], outlier_DT$fst), SNP=c(DT$id, paste0(outlier_DT$SNP,'_outlier'))) 
  dat$chr <- as.numeric(dat$chr)
  dat$mid_pos <- as.numeric(dat$mid_pos)
  dat$fst <- as.numeric(dat$fst)
  
  manhattan(dat, chr="chr",bp="mid_pos",p="fst", highlight1 = paste0(outlier_DT[outlier_DT$ed-outlier_DT$st==1000,4],'_outlier'), highlight2=paste0(outlier_DT[outlier_DT$ed-outlier_DT$st!=1000,4],'_outlier'), logp=FALSE, cex.axis = 1, ylim = c(0, 1.02), #subset(dat, chr == 6), xlim = c(21400000, 21420000),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="Fst", cex.lab=1.5, main = paste0(titlename," union outliers based 1000, 5000, 10000 bp window tests"), cex.main=1.5)
  dev.off()
  
  # output bed file
  colnames(outlier_DT)=c('Chromosome', 'Start', 'End', 'Feature', 'Fst')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    outlier_DT$Chromosome[outlier_DT$Chromosome==i] = chr_str_list[i]
  for(win in windows)
    outlier_DT$Feature[outlier_DT$End - outlier_DT$Start == win] = win
  write.table(outlier_DT, file = paste0(headname,"sliding.final.outlier.igv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}


############################## formal run ##############################
############################## formal run ##############################
############################## formal run ##############################

hudson_fst_sliding("SL_OBOYS2_noinvers.", "SL vs OBOYS2")
hudson_fst_sliding("SL_LOLA_noinvers.", "SL vs LOLA")
hudson_fst_sliding("NEH_UMFS_noinvers.", "NEH vs UMFS")
hudson_fst_sliding("CS_UMFS_noinvers.", "CS vs UMFS")
hudson_fst_sliding("CS_NEH_noinvers.", "CS vs NEH")
hudson_fst_sliding("CS_HC_noinvers.", "CS vs HC")
hudson_fst_sliding("CS_DEBY_noinvers.", "CS vs DEBY")
hudson_fst_sliding("CL_OBOYS2_noinvers.", "CL vs OBOYS2")