
#######################
#  Adjust p-value    #
#######################

# load the p-values (2032113 SNPs, from global shared 1x coverage data)

ps_adj <- function(pname, delta_name, outname1, outname2) {
  #pname = "./REF19_CHR19_NB_HC_out_all_fish.txt" 
  #delta_name = "./alt_REF19_CHR19_NB_HC.txt"
  dat1 = read.delim(pname, header = FALSE, sep='\t')
  dat2 = read.delim(delta_name, header = FALSE, sep='\t')
  # the code block tries to convert the number to chromosome IDs
  # chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  # for(i in seq(10)) 
  #   dat1$V1[dat1$V1==i] = chr_str_list[i]
  # for(i in seq(10)) 
  #   dat2$V1[dat2$V1==i] = chr_str_list[i]
  # dat1 <- dat1[with(dat1, order(V1, V2)),]
  # dat2 <- dat2[with(dat2, order(V1, V2)),]
  ps = dat1$V3
  id = paste0(dat1$V1,'_',dat1$V2)
  ps_cut = 0.05
  outlier1 <- as.data.frame(cbind(dat1$V1[ps<ps_cut], dat1$V2[ps<ps_cut], id[ps<ps_cut], dat1$V3[ps<ps_cut], dat2$V3[ps<ps_cut], dat2$V4[ps<ps_cut]))
  colnames(outlier1) <- c("chr","pos", "id", "p_value", "dp1", "dp2")
  outlier1 <- outlier1[with(outlier1, order(chr, pos)),]
  write.table(outlier1,outname1, row.names = FALSE, sep="\t", quote = FALSE)
  message("Number of outliers with p< 0.05:" ,length(outlier1$id))
  share1 = sum(sign(as.numeric(outlier1$dp1)) == sign(as.numeric(outlier1$dp2)))
  message("Number of outliers (p< ", ps_cut, ") with the same directionality: " ,share1)
  adj = p.adjust(ps, method = 'BH')
  #adj = p.adjust(ps, method = 'bonferroni')
  # check how many SNPs with FDR < 0.05
  adj_cut = 0.01
  length(adj[adj<adj_cut])
  outlier2 <- as.data.frame(cbind(dat1$V1[adj<adj_cut], dat1$V2[adj<adj_cut], id[adj<adj_cut], dat1$V3[adj<adj_cut], adj[adj<adj_cut], dat2$V3[adj<adj_cut], dat2$V4[adj<adj_cut]))
  colnames(outlier2) <- c("chr","pos", "id", "p_value","fdr", "dp1", "dp2")
  outlier2 <- outlier2[with(outlier2, order(chr, pos)),]
  write.table(outlier2,outname2, row.names = FALSE, sep="\t", quote = FALSE)
  message("Number of outliers after FDR:" ,length(outlier2$id))
  share2 = sum(sign(as.numeric(outlier2$dp1)) == sign(as.numeric(outlier2$dp2)))
  message("Number of outliers (p< ", adj_cut, ") with the same directionality: " ,share2)
}

ps_adj("REF19_CHR19_NB_HC_out_all_fish.txt", "alt_REF19_CHR19_NB_HC.txt","REF19_CHR19_NB_HC_ps_0.05.txt", "REF19_CHR19_NB_HC_adj_0.01.txt")
ps_adj("REF19_CHR19_SR_HC_out_all_fish.txt", "alt_REF19_CHR19_SR_HC.txt","REF19_CHR19_SR_HC_ps_0.05.txt", "REF19_CHR19_SR_HC_adj_0.01.txt")
ps_adj("REF19_SR_ARN_COH_out_all_fish.txt", "alt_REF19_SR_ARN_COH.txt","REF19_SR_ARN_COH_ps_0.05.txt", "REF19_SR_ARN_COH_adj_0.01.txt")


ps_adj("./Del20_data/ps_Del20_challenge.txt", "Del20_SGS_outlier.list")
ps_adj("./Del19_data/ps_Del19_challenge.txt", "Del19_SGS_outlier.list")
ps_adj("./WILD_HC_NB_data/ps_HC_NB_challenge.txt", "HC_NB_SGS_outlier.list")
ps_adj("./WILD_HC_SR_data/ps_HC_SR_challenge.txt", "HC_SR_SGS_outlier.list")


ps <- function(pname, outname) {
  #pname = "./REF19_SR_ARN_COH_out_all_fish.txt" 
  dat = read.delim(pname, header = FALSE, sep='\t')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat <- dat[with(dat, order(V1, V2)),]
  ps = dat$V3
  id = paste0(dat$V1,'_',dat$V2)
  hist(ps)
  # check how many SNPs with p-value < 0.05
  length(ps[ps<0.05])
  outlier <- as.data.frame(cbind(dat$V1[ps<0.05], dat$V2[ps<0.05], dat$V3[ps<0.05], id[ps<0.05]))
  colnames(outlier) <- c("chr","pos","p_value", "id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  write.table(outlier,outname, row.names = FALSE, sep="\t", quote = FALSE)
  message("Number of outliers:" ,length(outlier$id))
}

ps("./REF19_CHR19_NB_HC_out_all_fish.txt", "REF19_CHR19_NB_HC_fish_outlier.list")

ps("./REF19_CHR19_SR_HC_out_all_fish.txt", "REF19_CHR19_SR_HC_fish_outlier.list")
ps("./REF19_CHR19_REF20_CHR20_out_all_fish.txt", "REF19_CHR19_REF20_CHR20_fish_outlier.list")
ps("./REF20_CHR20_SR_HC_out_all_fish.txt", "REF20_CHR20_SR_HC_fish_outlier.list")

ps("./REF19_REF20_CHR19_CHR20_out_all_fish.txt", "REF19_REF20_CHR19_CHR20_fish_outlier.list")
ps("./REF19_SR_ARN_COH_out_all_fish.txt", "REF19_SR_ARN_COH_fish_outlier.list")


#######################
#  Shared outliers    #
#######################

make_ID <- function(file_name){
  dat = read.delim(file_name, header = TRUE, sep='\t')
  dat = dat[with(dat, order(chr, pos)),]
  dat_ID = dat$id
  message(length(dat$id))
  return(dat_ID)
}

HC_NB = make_ID('REF19_CHR19_NB_HC_fish_outlier.list')
Control = make_ID('REF19_SR_ARN_COH_fish_outlier.list')

REF19_CHR19_SR_HC = make_ID('REF19_CHR19_SR_HC_fish_outlier.list')
REF20_CHR20_SR_HC = make_ID('REF20_CHR20_SR_HC_fish_outlier.list')
REF19_CHR19_REF20_CHR20 = make_ID('REF19_CHR19_REF20_CHR20_fish_outlier.list')

a <- intersect(REF19_CHR19_SR_HC, REF20_CHR20_SR_HC)
length(a)

b <- intersect(REF19_CHR19_REF20_CHR20, intersect(REF19_CHR19_SR_HC,REF20_CHR20_SR_HC))
length(b)
write.table(b,"Fish_no_FDR_shared_3420.txt", row.names = FALSE, sep="\t", quote = FALSE)

dat = read.delim("GEA_outlier_226.txt", header = TRUE, sep='\t')
GEA_outlier_id = paste0(dat$chromo,'_',dat$position)

c <- intersect(b, GEA_outlier_id)
length(c)

############################################
#  check delta_p direction in outliers     #
############################################

setwd("~/Documents/HG/DelBay19_adult/08_indep_fish/shared_with_SGS")
# check outlier delta_p in each contrast
extract_outlier <- function(pname, snp_list){
  pname = "HC_NB_ps_outlier.list" # 2-side p-value
  dat = read.delim(pname, header = TRUE, sep='\t')
  # count how many SNPs with positive vs negative frequency change
  dat$id = paste0(dat$chr,'_',dat$pos)
  outlier <- dat[which(dat$id %in% snp_list),][,1:5] # modify the i to other common shared outliers
  colnames(outlier) <- c("chr","pos","p1","p2","deltap")
  cnt = length(outlier$deltap)
  pos_cnt = length(outlier[which(outlier$deltap>0),]$pos)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(outlier[which(outlier$deltap<0),]$pos)
  message("negative change SNPs: ", neg_cnt)
  #outlier1 <- outlier[with(outlier, order(-deltap)),] # order by deltap from high to low
  outlier1 <- outlier[with(outlier, order(chr, pos)),] # order by chromosome and position
  return(outlier1)
}

# count how many SNPs with negative frequency change
extract_neg <- function(pname){
  #pname = "indep_fish_NB_HC_ps_outlier.list" # 2-side p-value
  dat = read.delim(pname, header = TRUE, sep='\t')
  dat$id = paste0(dat$chr,'_',dat$pos)
  outlier <- dat[which(dat$id %in% i),][,c(1:5,8)] # modify the i to other common shared outliers
  colnames(outlier) <- c("chr","pos","p1","p2","deltap","id")
  cnt = length(outlier$deltap)
  message("total SNPs: ", cnt)
  neg_outlier = outlier[which(outlier$deltap<0),]
  neg_outlier_sort = neg_outlier[with(neg_outlier, order(chr, pos)),]
  message("negative change SNPs: ", length(neg_outlier_sort$chr))
  return(neg_outlier_sort)
}

# count how many SNPs with posative frequency change
extract_pos <- function(pname){
  #pname = "indep_fish_NB_HC_ps_outlier.list" # 2-side p-value
  dat = read.delim(pname, header = TRUE, sep='\t')
  dat$id = paste0(dat$chr,'_',dat$pos)
  outlier <- dat[which(dat$id %in% i),][,c(1:5,8)] # modify the i to other common shared outliers
  colnames(outlier) <- c("chr","pos","p1","p2","deltap","id")
  cnt = length(outlier$deltap)
  message("total SNPs: ", cnt)
  pos_outlier = outlier[which(outlier$deltap>=0),]
  pos_outlier_sort = pos_outlier[with(pos_outlier, order(chr, pos)),]
  message("posative change SNPs: ", length(pos_outlier_sort$chr))
  return(pos_outlier_sort)
}

Del19_fish_neg <- extract_neg("indep_fish_REF19_CHR19_ps_outlier.list")
Del19_fish_pos <- extract_pos("indep_fish_REF19_CHR19_ps_outlier.list")

HC_NB_fish_neg <- extract_neg("indep_fish_NB_HC_ps_outlier.list")
HC_NB_fish_pos <- extract_pos("indep_fish_NB_HC_ps_outlier.list")

HC_SR_fish_neg <- extract_neg("indep_fish_SR_HC_ps_outlier.list")
HC_SR_fish_pos <- extract_pos("indep_fish_SR_HC_ps_outlier.list")
