setwd("~/Documents/Ryan_workplace/DelBay_adult/08_fish_exact/Fish_exact_z_1x")
#######################
#  Adjust p-value    #
#######################

# load the p-values (2032113 SNPs, from global shared 1x coverage data)

ps_adj <- function(pname, outname) {
  #pname = "./REF19_SR_ARN_COH_out_all_fish.txt" 
  dat = read.delim(pname, header = FALSE, sep='\t')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat <- dat[with(dat, order(V1, V2)),]
  ps = dat$V3
  id = paste0(dat$V1,'_',dat$V2)
  hist(ps)
  adj = p.adjust(ps, method = 'BH')
  #adj = p.adjust(ps, method = 'bonferroni')
  hist(adj)
  # check how many SNPs with FDR < 0.05
  length(adj[adj<0.05])
  outlier <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], adj[adj<0.05], id[adj<0.05]))
  colnames(outlier) <- c("chr","pos","p_value","fdr", "id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  write.table(outlier,outname, row.names = FALSE, sep="\t", quote = FALSE)
  message("Number of outliers:" ,length(outlier$id))
}

ps_adj("REF19_CHR19_NB_HC_out_all_fish.txt", "REF19_CHR19_NB_HC_out_0.05_fish.txt")
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


