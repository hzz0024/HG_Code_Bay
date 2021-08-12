setwd("~/Documents/Ryan_workplace/CVreseq_hudson_fst_nobiginvers")
library(hash)
source("manhattan.R")
library(hash)
library(dgconstraint)
options(scipen=999)
library("gtools")

# headname = "CL_OBOYS2_noinvers."
# titlename = 'CL vs OBOYS2'
hudson_fst <- function(headname, titlename){
  #jpeg(paste0('C_score/',headname,"sliding.zfst.hudson.jpg"), width = 16, height = 9, units = 'in', res = 300)
  #jpeg(paste0('C_score/',headname,"zfst.hudson.jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(1,1))
  for(win in c(1000)){
    name = paste0(headname, win, "bp.", "s", win/5, ".csv")
    #name = paste0(headname, win, "bp.","csv")
    DT = read.delim(name, header = TRUE, sep=',')
    mid_pos <- round((DT$start + DT$end)/2)
    id = paste0(DT$scaffold,'_',mid_pos)
    DT <- as.data.frame(cbind(DT,mid_pos, id))
    #DT <- DT[complete.cases(DT), ]
    DT[,9][DT[,9]<0] = 0.000001 #@chnage
    zfst <- (DT[,9] - mean(DT[,9], na.rm=T))/(sd(DT[,9], na.rm = T))# @change
    dat <- data.frame(chr=DT$scaffold, start=DT$start, end=DT$end, mid_pos=DT$mid_pos, SNP=DT$id, fst=DT[,9], zfst = zfst) # @change
    dat$chr <- as.numeric(dat$chr)
    dat$mid_pos <- as.numeric(dat$mid_pos)
    dat$fst <- as.numeric(dat$fst)
    dat$zfst <- as.numeric(dat$zfst)
    # calculate the 99.9% quantile
    dat1 <- dat[complete.cases(dat), ]
    thred=quantile(dat1$zfst, 0.999, na.rm=T)
    cnt = length(dat1[dat1$zfst>thred[[1]],1])
    outlier = dat1[dat1$zfst>thred[[1]],5]
    #cnt = sum(dat$zfst>5)
    #outlier = dat[dat$zfst>5,5]
    print(paste0("Number of outlier in ", titlename," at 99.9% quantile is ", cnt, " (window size ", win, ")" ))
    #manhattan(dat1, chr="chr",bp="mid_pos",p="zfst", highlight1 = outlier, logp=FALSE, cex.axis = 1, ylim = c(0, max(dat$zfst, na.rm=T)+0.2), #subset(dat, chr == 8)
    #          col=c("grey","black"),genomewideline=F, suggestiveline=F,
    #          ylab="ZFst", cex.lab=1.5, main = paste0(titlename, " ZFst window size:", win, ' bp'), cex.main=1.5)
    #dev.off()
    sliding_DT = data.frame(Chromosome=dat$chr, Start=dat$start, End=dat$end, Fst=dat$fst, ZFst=dat$zfst, outlier=dat$zfst>thred[[1]]) #@changed at 12/10/2020 to output the ZFst for each window
    write.table(sliding_DT, file = paste0('C_score/',headname,"all.",win, "bp.txt"), sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE) #@changed at 12/10/2020 to output the ZFst for each window
    }
}

hudson_fst("SL_OBOYS2_noinvers.", "SL vs OBOYS2")
hudson_fst("SL_LOLA_noinvers.", "SL vs LOLA")
hudson_fst("NEH_UMFS_noinvers.", "NEH vs UMFS")
hudson_fst("CS_UMFS_noinvers.", "CS vs UMFS")
hudson_fst("CS_NEH_noinvers.", "CS vs NEH")
hudson_fst("CS_DEBY_noinvers.", "CS vs DEBY")
hudson_fst("CL_OBOYS2_noinvers.", "CL vs OBOYS2")

hudson_fst("CS_HC_noinvers.", "CS vs HC")
hudson_fst("HCVA_CLP_noinvers.", "HCVA vs CLP")
hudson_fst("HC_CLP_noinvers.", "HC vs CLP")
hudson_fst("CS_HCVA_noinvers.", "CS vs HCVA")

win = 1000
heads = c("SL_OBOYS2_noinvers.","SL_LOLA_noinvers.","NEH_UMFS_noinvers.","CS_UMFS_noinvers.",'CS_DEBY_noinvers.', 'CS_NEH_noinvers.','CL_OBOYS2_noinvers.',"CS_HC_noinvers.")
combs = combinations(8, 2, 1:8)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0('C_score/',headname1,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0('C_score/',headname2,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  
  DT1$ZFst[which(DT1$outlier==FALSE)] = 0
  DT2$ZFst[which(DT2$outlier==FALSE)] = 0
  
  DT1$ZFst[which(DT1$outlier==TRUE)] = 1
  DT2$ZFst[which(DT2$outlier==TRUE)] = 1

  del_idx = c(which(is.na(DT1$ZFst)), which(is.na(DT2$ZFst)))
  
  re1 = single_c_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  print(paste0("C-score is ",re1))
  re2 = single_p_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  print(paste0("p-value is ",re2))
  re3 = single_c_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("C-score is ",re3))
  re4 = single_p_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("p-value is ",re4))
}

########################### C score running for target contrasts ##################################

win = 1000
heads = c("CS_HC_noinvers.",'HCVA_CLP_noinvers.', 'HC_CLP_noinvers.','CS_HCVA_noinvers.')
combs = combinations(4, 2, 1:4)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0('C_score/',headname1,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0('C_score/',headname2,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  
  DT1$ZFst[which(DT1$outlier==FALSE)] = 0
  DT2$ZFst[which(DT2$outlier==FALSE)] = 0
  
  DT1$ZFst[which(DT1$outlier==TRUE)] = 1
  DT2$ZFst[which(DT2$outlier==TRUE)] = 1
  
  del_idx = c(which(is.na(DT1$ZFst)), which(is.na(DT2$ZFst)))
  
  re1 = single_c_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  print(paste0("C-score is ",re1))
  re2 = single_p_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  print(paste0("p-value is ",re2))
  re3 = single_c_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("C-score is ",re3))
  re4 = single_p_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("p-value is ",re4))
}

######################## all contrasts ######################## 
win = 1000
heads = c("CS_UMFS_noinvers.",'CS_DEBY_noinvers.', 'CS_NEH_noinvers.')
combs = combinations(3, 2, 1:3)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0('C_score/',headname1,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0('C_score/',headname2,"all.",win, "bp.txt"), header = TRUE, sep='\t')
  
  DT1$ZFst[which(DT1$outlier==FALSE)] = 0
  DT2$ZFst[which(DT2$outlier==FALSE)] = 0
  
  DT1$ZFst[which(DT1$outlier==TRUE)] = 1
  DT2$ZFst[which(DT2$outlier==TRUE)] = 1
  
  del_idx = c(which(is.na(DT1$ZFst)), which(is.na(DT2$ZFst)))
  
  #re1 = single_c_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  #print(paste0("C-score is ",re1))
  #re2 = single_p_chisq(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], num_permute=10000, na.rm = F)
  #print(paste0("p-value is ",re2))
  re3 = single_c_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("C-score is ",re3))
  re4 = single_p_hyper(DT1$ZFst[-del_idx], DT2$ZFst[-del_idx], na.rm = F)
  print(paste0("p-value is ",re4))
}
