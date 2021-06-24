setwd("~/Documents/Ryan_workplace/DelBay_adult/11_SGS")
#######################
#  Adjust p-value    #
#######################

# load the p-values (2032113 SNPs, from global shared 1x coverage data)

ps_adj <- function(pname, outname) {
  #pname = "barnard_REF19_CHR19.txt" # 2-side p-value
  dat = read.delim(pname, header = FALSE, sep='\t')
  chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
  for(i in seq(10)) 
    dat$V1[dat$V1==i] = chr_str_list[i]
  dat[with(dat, order(V1, V2)),]
  ps = dat$V6
  ps_cnt = length(ps[ps<0.05])
  message("p-value < 0.05:", ps_cnt)
  id = paste0(dat$V1,'_',dat$V2)
  
  hist(ps)
  adj = p.adjust(ps, method = 'BH')
  hist(adj)
  # check how many SNPs with FDR < 0.05
  length(adj[adj<0.05])
  outlier <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], adj[adj<0.05], id[adj<0.05]))
  colnames(outlier) <- c("chr","pos","p_value","fdr", "id")
  outlier <- outlier[with(outlier, order(chr, pos)),]
  write.table(outlier,outname, row.names = FALSE, sep="\t", quote = FALSE)
  message("FDR < 0.05:" ,length(outlier$id))
}

ps_adj("barnard_REF19_CHR19.txt", "Del19_outlier.list")
ps_adj("barnard_REF20_CHR20.txt", "Del20_outlier.list")
ps_adj("barnard_SR_HC.txt", "SR_HC_outlier.list")
ps_adj("barnard_NB_HC.txt", "NB_HC_outlier.list")

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

Del19 = make_ID('Del19_outlier.list')
Del20 = make_ID('Del20_outlier.list')
HC_SR = make_ID('SR_HC_outlier.list')
HC_NB = make_ID('NB_HC_outlier.list')

a <- intersect(Del19, Del20)
length(a)
b <- intersect(Del19, HC_SR)
length(b)
c <- intersect(Del19, HC_NB)
length(c)
e <- intersect(Del20, HC_SR)
length(e)
d <- intersect(Del20, HC_NB)
length(d)
f <- intersect(HC_NB, HC_SR)
length(f)

