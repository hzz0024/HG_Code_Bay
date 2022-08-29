
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/11_SGS")

name1 = "REF19_CHR19_REF20_CHR20_FDR0.05_same.txt"
DT1 = read.delim(name1, header = TRUE, sep='\t')
DT1$p_c[DT1$p_c == 0] = 0.00001
DT1$p_w[DT1$p_w == 0] = 0.00001
id1 = paste0(DT1$chromo,'_',DT1$position)
message("Number of oultier SNPs in ", name1, " is ", length(id1))
name2 = "REF19_CHR19_NB_HC_FDR0.05_same.txt"
DT2 = read.delim(name2, header = TRUE, sep='\t')
DT2$p_c[DT2$p_c == 0] = 0.00001
DT2$p_w[DT2$p_w == 0] = 0.00001
id2 = paste0(DT2$chromo,'_',DT2$position)
message("Number of oultier SNPs in ", name2, " is ", length(id2))
idx1 = (id1 %in% id2)
message("Number of shared SNPs is ", sum(idx1))
shared = DT1[idx1,]

write.table(shared, "./shared_snps.txt", row.names=F, quote=F, sep="\t")

idx2 = (id2 %in% id1)
shared = DT2[idx2,]
message("Number of shared SNPs is ", sum(idx2))

bed <- data.frame(shared[1], shared[2], shared[2])
bed <- data.frame(shared[1], shared[2]-1000, shared[2]+1000)
write.table(bed, "./shared_bed_188_2K.txt", row.names=F, quote=F, sep="\t")

#########################################
## shared SNPs test for directionality ## 
#########################################  
check_direct <- function(dat1){
  #pname = "./REF19_CHR19_NB_HC_out_all_fish.txt"
  #dat = read.delim(pname, header = TRUE, sep='\t')
  dp1 = dat1$dp_c
  dp2 = dat1$dp_w
  share_all = sum(sign(as.numeric(dp1)) == sign(as.numeric(dp2)))
  message("number of SNPs with same directionality is ", share_all, " out of ", dim(dat1)[1], " (", share_all/dim(dat1)[1], ")")
}

check_direct(shared_challenge, shared_wild)
check_direct(shared_wild)

#########################################
## random SNPs for ps                  ## 
######################################### 
pname = "./ps_Del19_challenge.txt"
pname = "./ps_Del20_challenge.txt"
pname = "./ps_Del19_HC_NB.txt"
dat = read.delim(pname, header = FALSE, sep='\t')
head(dat)
random <- sample(dat$V4,5000,replace=FALSE)
write.table(random, "./random_ref19.txt", row.names=F, quote=F, sep="\t")
write.table(random, "./random_ref20.txt", row.names=F, quote=F, sep="\t")
write.table(random, "./random_NB.txt", row.names=F, quote=F, sep="\t")
