setwd("./")
pname = "ps_Del19_challenge.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
message("number of SNPs with positive delta_p is ", length(dat$V5[which(dat$V5>=0)]))
dat$delta_p <- dat$V3-dat$V4
dat$SNP = paste0(dat$V1,'_',dat$V2)
dat$V6[dat$V6 == 0] = 0.00001
dat$adj = p.adjust(dat$V6, method = 'BH')
colnames(dat)=c('chromo', 'position', 'p1', 'p0', 'D', 'ps', 'raw_candidates', 'delta_p', 'SNP', 'adj')
dim(dat)

idx = dat$position>0 # all TRUE
for(i in seq(dim(dat)[1])){
#for(i in seq(1:10000)){
  chr = dat$chromo[i]
  pos = dat$position[i]
  cnt = 0
  # pre range
  for(j in seq(i,i-10)){
    if(j<1||j>dim(dat)[1])
      break
    cur_chr = dat$chromo[j]
    cur_pos = dat$position[j]
    if(cur_chr==chr && cur_pos>(pos-250))
      cnt = cnt + 1
  }
  # post range
  for(j in seq(i,i+10)){
    if(j<1||j>dim(dat)[1])
      break
    cur_chr = dat$chromo[j]
    cur_pos = dat$position[j]
    if(cur_chr==chr && cur_pos<(pos+250))
      cnt = cnt + 1
  }
  if(cnt<13)
    idx[i] = FALSE
}

dat1 =dat[idx,]

