
find_shared <- function(f1_name,f2_name) {
  #f1_name = "ps_Del19_challenge.txt"
  file1 = f1_name
  df1 <- read.delim(file1, header = FALSE, sep='\t')
  out1 = df1[which(df1$V3 < 0.05),]
  out1$id = paste0(out1$V1,'_',out1$V2)
  message("File ", f1_name, " outlier: ", length(out1$id))
  #f2_name = "barnard_REF19_CHR19.txt"
  file2 = f2_name
  df2 <- read.delim(file2, header = FALSE, sep='\t')
  out2 = df2[which(df2$V6 < 0.05),]
  out2$id = paste0(out2$V1,'_',out2$V2)
  message("File ", f2_name, " outlier: ", length(out2$id))
  #Check shared SNPs
  share_snps <- out2[(out2$id %in% out1$id),][,1:5]
  message("Shared SNPs: ", length(share_snps$V1)) 
  pos_cnt = length(share_snps[which(share_snps$V5>0),]$V1)
  message("positive change SNPs: ", pos_cnt)
  neg_cnt = length(share_snps[which(share_snps$V5<0),]$V1)
  message("negative change SNPs: ", neg_cnt)
  return(share_snps)
 }

Del19 <- find_shared("ps_Del19_challenge.txt", "barnard_REF19_CHR19.txt")
Del20 <- find_shared("ps_Del20_challenge.txt", "barnard_REF20_CHR20.txt")
HC_NB <- find_shared("ps_HC_NB_challenge.txt", "barnard_NB_HC.txt")
HC_NB <- find_shared("ps_HC_SR_challenge.txt", "barnard_SR_HC.txt")

n_bootstraps=1000
boot_cnt = rep(0, n_bootstraps)
a = seq(1,2032113,1)
for (i in 1:n_bootstraps){
  o1 = sample(a, 3265)
  o2 = sample(a, 2559)
  cnt = length(intersect(o1,o2))
  boot_cnt[i] = cnt
}

b=t.test(boot_cnt)

hist(boot_cnt)
