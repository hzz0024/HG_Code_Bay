
#######################
#  Adjust p-value    #
#######################

# load the p-values (2032113 SNPs, from global shared 1x coverage data)

ps_adj <- function(pname, outname) {
  #pname = "ps_Del20_challenge.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  ps = dat$V3
  id = paste0(dat$V1,'_',dat$V2)
  hist(ps)
  adj = p.adjust(ps, method = 'BH')
  hist(adj)
  # check how many SNPs with FDR < 0.05
  length(adj[adj<0.05])
  outlier <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], adj[adj<0.05], id[adj<0.05]))
  colnames(outlier) <- c("chr","pos","p_value","fdr", "id")
  outlier[with(outlier, order(chr, pos)),]
  write.table(outlier,outname, row.names = FALSE, sep="\t", quote = FALSE)
  message("Number of outliers:" ,length(outlier$id))
}

ps_adj("ps_Del20_challenge.txt", "Del20_SGS_outlier.list")
ps_adj("ps_Del19_challenge.txt", "Del19_SGS_outlier.list")
#######################
#  Extract outliers   #
#######################

# extract the outlier for snp list creation
fname1 = 'Del19_20_1x_share_snps.list'  
dat_sites = read.delim(fname1, header = FALSE, sep='\t')
dat_sites$id = paste0(dat_sites$V1,'_',dat_sites$V2)
colnames(dat_sites) <- c("chr","pos","major","minor", "id")
outlier_list <- as.data.frame(cbind(dat_sites$chr[dat_sites$id %in% outlier$id], 
                                    dat_sites$pos[dat_sites$id %in% outlier$id], 
                                    dat_sites$major[dat_sites$id %in% outlier$id], 
                                    dat_sites$minor[dat_sites$id %in% outlier$id]))
write.table(outlier_list,"Del19_SGS_outlier_3006.list", row.names = FALSE, sep="\t", quote = FALSE)

#######################
#  Shared outliers    #
#######################

make_ID <- function(file_name){
  dat = read.delim(file_name, header = FALSE, sep='\t')
  dat = dat[with(dat, order(V1, V2)),]
  adj = p.adjust(dat$V3, method = 'BH')
  datt = dat
  datt = cbind(dat, adj)
  idx = datt$adj < 0.01
  message(paste0('0.01: ',length(datt$adj[idx])))
  dat_ID = paste0(datt$V1[idx],'_',datt$V2[idx])
  return(dat_ID)
}

REF_SR = make_ID('SGS_REF_SR_2side_ps.txt')
HC_NB = make_ID('SGS_HC_NB_2side_ps.txt')
HC_SR = make_ID('SGS_HC_SR_2side_ps.txt')
CH_REF = make_ID('SGS_CH_REF_2side_ps.txt')

a <- intersect(CH_REF, HC_SR)
length(a)

b <- intersect(CH_REF, HC_NB)
length(b)

c <- intersect(HC_NB, HC_SR)
length(c)

d <- intersect(HC_NB, REF_SR)
length(d)

e <- intersect(HC_SR, REF_SR)
length(e)

f <- intersect(CH_REF, REF_SR)
length(f)
