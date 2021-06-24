
find_shared <- function(f1_name,f2_name) {
  #f1_name = "Del19_SGS_outlier.list"
  file1 = f1_name
  df1 <- read.delim(file1, header = TRUE, sep='\t')
  message("File ", f1_name, " SNPs: ", length(df1$id))
  #f2_name = "HC_NB_SGS_outlier.list"
  file2 = f2_name
  df2 <- read.delim(file2, header = TRUE, sep='\t')
  message("File ", f2_name, " SNPs: ", length(df2$id))
  share_snps <- df1[(df1$id %in% df2$id),][,1:5]
  message("Shared SNPs: ", length(share_snps$id))
  return(share_snps)
 }

find_shared("Del19_SGS_outlier.list", "HC_NB_SGS_outlier.list")
find_shared("Del19_SGS_outlier.list", "HC_SR_SGS_outlier.list")

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
