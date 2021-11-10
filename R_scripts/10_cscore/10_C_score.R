################################################
###   format the indepedent SGS results      ###
################################################

setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/10_score")
library(hash)
library(dgconstraint)
options(scipen=999)
library("gtools")

ps_adj05 <- function(pname, outname2) {
  # load the p-values (2335739 SNPs, from global snp list
  #pname = "ps_Del19_challenge.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  # replace chromosome if it is numerical
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  # count how many SNPs with p < 0.05
  ps = dat$V6
  ps[which(ps==0)]=0.00001
  ps_cnt = length(ps[ps<0.05])
  message("p-value < 0.05:", ps_cnt)
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$V1,'_',dat$V2)
  ps_outlier <- dat[which(dat$V6<0.05),][,1:5]
  pos_cnt = length(ps_outlier[which(ps_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(ps_outlier[which(ps_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  # adjust the p-values with FDR correction
  dat$adj = p.adjust(ps, method = 'BH')
  dat$bi[which(dat$adj>=0.05)] = 0
  dat$bi[which(dat$adj<0.05)] = 1
  # count how many SNPs with FDR < 0.05
  message("FDR < 0.05:" ,length(outlier2$id))
  # count how many SNPs with positive vs negative frequency change
  adj_outlier <- dat[which(dat$adj<0.05),][,1:5]
  pos_cnt = length(adj_outlier[which(adj_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(adj_outlier[which(adj_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  
  # output a table for SNPs with FDR < 0.05
  outlier2 <- as.data.frame(cbind(dat$V1, dat$V2, dat$adj, dat$bi))
  colnames(outlier2) <- c("chr","pos","FDR","bi")
  outlier2 <- outlier2[with(outlier2, order(chr, pos)),]
  write.table(outlier2,outname2, row.names = FALSE, sep="\t", quote = FALSE)
}


ps_adj01 <- function(pname, outname2) {
  # load the p-values (2335739 SNPs, from global snp list
  #pname = "ps_Del19_challenge.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  # replace chromosome if it is numerical
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  # count how many SNPs with p < 0.01
  ps = dat$V6
  ps[which(ps==0)]=0.00001
  ps_cnt = length(ps[ps<0.01])
  message("p-value < 0.01:", ps_cnt)
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$V1,'_',dat$V2)
  ps_outlier <- dat[which(dat$V6<0.01),][,1:5]
  pos_cnt = length(ps_outlier[which(ps_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(ps_outlier[which(ps_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  # adjust the p-values with FDR correction
  dat$adj = p.adjust(ps, method = 'BH')
  dat$bi[which(dat$adj>=0.01)] = 0
  dat$bi[which(dat$adj<0.01)] = 1
  # count how many SNPs with FDR < 0.01
  message("FDR < 0.01:" ,length(outlier2$id))
  # count how many SNPs with positive vs negative frequency change
  adj_outlier <- dat[which(dat$adj<0.01),][,1:5]
  pos_cnt = length(adj_outlier[which(adj_outlier$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(adj_outlier[which(adj_outlier$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  
  # output a table for SNPs with FDR < 0.01
  outlier2 <- as.data.frame(cbind(dat$V1, dat$V2, dat$adj, dat$bi))
  colnames(outlier2) <- c("chr","pos","FDR","bi")
  outlier2 <- outlier2[with(outlier2, order(chr, pos)),]
  write.table(outlier2,outname2, row.names = FALSE, sep="\t", quote = FALSE)
}

ps_adj05("ps_Del19_challenge.txt", "REF19_CHR19_FDR_outlier.list")
ps_adj05("ps_Del19_HC_NB.txt", "NB_HC_FDR_outlier.list")
ps_adj05("ps_Del19_HC_SR.txt", "SR_HC_FDR_outlier.list")

ps_adj01("ps_Del19_challenge.txt", "REF19_CHR19_FDR_outlier.list")
ps_adj01("ps_Del19_HC_NB.txt", "NB_HC_FDR_outlier.list")
ps_adj01("ps_Del19_HC_SR.txt", "SR_HC_FDR_outlier.list")

###################################
###   for pairwise Cscore tests ###
###################################

heads = c("REF19_CHR19",'NB_HC', 'SR_HC')
combs = combinations(3, 2, 1:3)
for(i in 1:length(combs[,1])){
  print(combs[i,])
  print(heads[combs[i,1]])
  print(heads[combs[i,2]])
  headname1 = heads[combs[i,1]]
  DT1 = read.delim(paste0(headname1,"_FDR_outlier.list"), header = TRUE, sep='\t')
  headname2 = heads[combs[i,2]]
  DT2 = read.delim(paste0(headname2,"_FDR_outlier.list"), header = TRUE, sep='\t')
  #del_idx = c(which(is.na(DT1$bi)), which(is.na(DT2$bi)))
  DT = cbind(DT1$bi, DT2$bi)
  #re1 = single_c_chisq(DT1$bi[-del_idx], DT2$bi[-del_idx], num_permute=10000, na.rm = F)
  #print(paste0("C-score is ",re1))
  #re2 = single_p_chisq(DT1$bi[-del_idx], DT2$bi[-del_idx], num_permute=10000, na.rm = F)
  #print(paste0("p-value is ",re2))
  re3 = single_c_hyper(as.numeric(DT1$bi), as.numeric(DT2$bi), na.rm = T)
  print(paste0("C-score is ",re3))
  re4 = single_p_hyper(as.numeric(DT1$bi), as.numeric(DT2$bi), na.rm = F)
  print(paste0("p-value is ",re4))
  res5=prop_adaptive <- estimate_pa(DT, ndigits = 4, show.plot = F, na.rm = F)
  print(paste0("proportion of genes that can potentially contribute to adaptation ",res5))
}


max_likelihood <- likelihood_pa(prop_adaptive, copper, na.rm = F)

