#install_github('tavareshugo/windowscanr')
library(windowscanr)
source("manhattan.R")
####################################
##########  read input     ########
####################################
setwd("~/Dropbox/Mac/Documents/HG/DelBay19_adult/11_SGS/11_SGS_boot_window")
pname = "ps_Del19_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t',col.names=c('chromo', 'position', 'p1', 'p0', 'delta_p', 'ps', 'raw_candidates'))
dat = dat[order(dat[,1], dat[,2]),]
dat$ps[dat$ps == 0] = 0.00001
dat$adj = p.adjust(dat$ps, method = 'BH')
dat$SNP = paste0(dat$chromo,'_',dat$position)

# replace chromosome if it is numerical
chr_num_list = seq(10)
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(1:10)) 
  dat$chromo[dat$chromo==chr_str_list[i]] = chr_num_list[i]
dat$chromo = as.numeric(dat$chromo)

manhattan(chr="chromo",bp="position",p="adj", snp = "SNP", dat, logp=TRUE, cex.axis = 1.2,
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="2019 Surv - Ref SGS -log (FDR)", cex.lab=1) 
abline(h=-log(0.05, 10), lty = 2, col = "green", cex.lab=1.5)
abline(h=-log(0.01, 10), lty = 2, col = "red", cex.lab=1.5)
#idx = sort(sample(seq(1:dim(dat)[1]), 100))
#dat = dat[idx,]

mean_na_rm <- function(x, ...){
  mean(x, na.rm = TRUE, ...)
}

pos_win <- winScan(x = dat, 
                   groups = "chromo", 
                   position = "position", 
                   values = c("ps"), 
                   win_size = 500,
                   win_step = 250,
                   funs = c("mean_na_rm"))

save(pos_win, file = "Del19_challenge_pos_win.RData")

pos_win_noNA <- subset(pos_win, !is.na(ps_mean_na_rm))
pos_win_noNA
plot(-log(pos_win_noNA$ps_mean_na_rm, 10))
set.seed(1)
##step1. get the null distributions for counts of markers below a p theshold
##note that the p-values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values
get.bootstrap = function(DF, id, reps, cut){
  chromo = c()
  win_start = c()
  win_end = c()
  win_mid = c()
  ps_n = c()
  ps_mean_na_rm = c()
  ps = c()
  for (i in 1:dim(DF)[1]){
    chr_ = DF$chromo[i]
    ws = DF$win_start[i]
    we = DF$win_end[i]
    wm = DF$win_mid[i]
    snp_cnt = DF$ps_n[i]
    ps_mean = DF$ps_mean_na_rm[i]
    sub = sample(DF$ps_mean_na_rm, reps)
    obs = DF$ps_mean_na_rm[i]
    extreme.count = length(sub[sub < obs])
    p.value = extreme.count/reps
    chromo = c(chromo, chr_)
    win_start = c(win_start, ws)
    win_end = c(win_end, we)
    win_mid = c(win_mid, wm)
    ps_n = c(ps_n, snp_cnt)
    ps_mean_na_rm = c(ps_mean_na_rm, ps_mean)
    ps = c(ps, p.value)
    # if(p.value < cut){
    #   chromo = c(chromo, chr_)
    #   win_start = c(win_start, ws)
    #   win_end = c(win_end, we)
    #   win_mid = c(win_mid, wm)
    #   ps_n = c(ps_n, snp_cnt)
    #   ps_mean_na_rm = c(ps_mean_na_rm, ps_mean)
    #   ps = c(ps, p.value)
    # }
  }
  results = data.frame(chromo,win_start, win_end, win_mid, ps_n, ps_mean_na_rm, ps)
  #write.table(results, paste(GLOBAL_NAME, id, reps, cut, sep = "_"), row.names = F, quote = F)
  return(results)
}##Arugments: DF = the dataframe; id = a string to id the file from when exported


##step1. get the null distributions for counts of markers below a p theshold
##note that the p-values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values

set.seed(1)
dat <- get.bootstrap(pos_win_noNA, "test", 1000, 0.05)
write.table(dat, "Del19_challenge_pos_win_bootstrap.txt", row.names = F, col.names = T, quote = F)

name = "Del19_challenge_pos_win_bootstrap.txt"
dat = read.delim(name, header = T, sep=' ')
dat$ps[dat$ps == 0] = 0.001
dat = dat[which(dat$ps_n >=8),] # exclude windows with SNP counts < 8
dat$adj = p.adjust(dat$ps, method = "BH")
dat$SNP = paste0(dat$chromo,'_',dat$win_mid)

# replace chromosome if it is numerical
chr_num_list = seq(10)
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(1:10)) 
  dat$chromo[dat$chromo==chr_str_list[i]] = chr_num_list[i]
dat$chromo = as.numeric(dat$chromo)

manhattan(chr="chromo",bp="win_mid",p="ps", snp = "SNP", dat, logp=TRUE, cex.axis = 1.2,
          col=c("grey50","black"),genomewideline=F, suggestiveline=F,
          ylab="2019 Surv - Ref 500bp/win -log (p-value)", cex.lab=1) 
abline(h=-log(0.05, 10), lty = 2, col = "green", cex.lab=1.5)
abline(h=-log(0.01, 10), lty = 2, col = "red", cex.lab=1.5)

# check ps < 0.01
dat1 = dat[which(dat$ps < 0.01),]
message(paste0("total number of outlier window is ", dim(dat1)[1]))
dat_ = dat1[with(dat1, order(chromo, win_mid)),]
bed_list <- paste0(dat_$chromo, "\t" , dat_$win_start, "\t", dat_$win_end)
print(head(bed_list))
write.table(bed_list, paste0("CHR19_REF19_win_p01.bed"), row.names=F, col.names = F, quote=F, sep="\t")






dat$SNP = paste0(dat$chromo,'_',dat$win_mid)
dat1$SNP = paste0(dat1$chromo,'_',dat1$win_mid)
length(intersect(dat$SNP, dat1$SNP))

# replace chromosome if it is numerical
chr_num_list = seq(10)
chr_str_list = c('NC_035780.1','NC_035781.1','NC_035782.1','NC_035783.1','NC_035784.1','NC_035785.1','NC_035786.1','NC_035787.1','NC_035788.1','NC_035789.1')
for(i in seq(1:10)) 
  pos_win_noNA$chromo[pos_win_noNA$chromo==chr_str_list[i]] = chr_num_list[i]
pos_win_noNA$chromo = as.numeric(pos_win_noNA$chromo)
pos_win_noNA$SNP = paste0(pos_win_noNA$chromo,'_',pos_win_noNA$position)
is.odd <- function(x) x %% 2 != 0
colors = is.odd(pos_win_noNA$chromo)

colors[colors == TRUE] <- 'black'
colors[colors == FALSE] <- 'grey'

#make significant points red
threshold = 0.05
x = dat$adj < threshold
for (i in 1:length(x)){
  if (x[i] == TRUE){
    colors[i] <- 'red'
  }
}





plot(-log(ps_mean_na_rm, 10)~win_mid, data = pos_win_noNA, main = NULL, pch = 19, col = colors, cex = 0.2, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(fdr)'))
axis(1, at = 1000000, labels = F)
axis(2, labels = NULL, las = 1, cex.axis = 1)

axis(2, labels = NULL, las = 1, cex.axis = 1)

mtext(lg.names, side = 1, at = mids, line = .75, cex = 1)

segments(0, -0.1, 100000, -0.1, lty = 1, col = "red", lwd = 8)
