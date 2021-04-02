#######################
#  Adjust p-value    #
#######################
# load the p-values (~1.7 million SNPs, old dataset)
fname = 'ps_1side.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat[with(dat, order(V1, V2)),]
ps = dat$V3
# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.05], dat$V2[adj<0.05], dat$V3[adj<0.05], adj[adj<0.05]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"ps_1side_adj.txt", row.names = FALSE, sep="\t", quote = FALSE)

#######################################
#  Adjust p-value for 2-side tests    #
#######################################
# load the p-values (~1.7 million SNPs, old dataset)
fname = 'ps_2side.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat[with(dat, order(V1, V2)),]
ps = dat$V3

# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.01], dat$V2[adj<0.01], dat$V3[adj<0.01], adj[adj<0.01]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"ps_2side_adj.txt", row.names = FALSE, sep="\t", quote = FALSE)



#######################################
#  Adjust p-value for CH_REF tests    #
#######################################
fname = 'SGS_CH_REF_2side_ps.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
ps = dat$V3

# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.01], dat$V2[adj<0.01], dat$V3[adj<0.01], adj[adj<0.01]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"SGS_CH_REF_fdr0.01.txt", row.names = FALSE, sep="\t", quote = FALSE)

#######################################
#  Adjust p-value for HC_NB tests    #
#######################################
fname = 'SGS_HC_NB_2side_ps.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
ps = dat$V3

# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.01], dat$V2[adj<0.01], dat$V3[adj<0.01], adj[adj<0.01]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"SGS_HC_NB_fdr0.01.txt", row.names = FALSE, sep="\t", quote = FALSE)


#######################################
#  Adjust p-value for HC_SR tests    #
#######################################
fname = 'SGS_HC_SR_2side_ps.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
ps = dat$V3

# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.01], dat$V2[adj<0.01], dat$V3[adj<0.01], adj[adj<0.01]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"SGS_HC_SR_fdr0.01.txt", row.names = FALSE, sep="\t", quote = FALSE)


#######################################
#  Adjust p-value for REF_SR tests    #
#######################################
fname = 'SGS_REF_SR_2side_ps.txt'  
dat = read.delim(fname, header = FALSE, sep='\t')
dat = dat[with(dat, order(V1, V2)),]
ps = dat$V3

# load SNP dataset
library(qvalue)
# check how many SNPs with p_value < 0.05
length(which(ps<0.05))
# plot the p-value histogram
hist(ps)
# start p_value adjustment
adj = p.adjust(ps, method = 'BH')
# check how many SNPs with FDR < 0.2
length(adj[adj<0.01])
outlier <- as.data.frame(cbind(dat$V1[adj<0.01], dat$V2[adj<0.01], dat$V3[adj<0.01], adj[adj<0.01]))
colnames(outlier) <- c("chr","pos","p_value","fdr")
outlier
# export the dataframe to txt file
write.table(outlier,"SGS_REF_SR_fdr0.01.txt", row.names = FALSE, sep="\t", quote = FALSE)

################# check if any shared outliers #################

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
