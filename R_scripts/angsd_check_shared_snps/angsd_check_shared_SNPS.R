setwd("~/Documents/Ryan_workplace/DelBay_adult/batch_effect/downsampling")
# this script is used to find common shared snps between two snp list files
file1 = "Del19_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt"
df1 <- read.delim(file1, header = FALSE, sep='\t')
df1$V5 <- paste(df1$V1, df1$V2, sep="_")
length(df1$V5)
#file2 = "Del20_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt"
file2 = "down_sampling_0p7x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt"
df2 <- read.delim(file2, header = FALSE, sep='\t')
df2$V5 <- paste(df2$V1, df2$V2, sep="_")
length(df2$V5)
share_snps <- df1[(df1$V5 %in% df2$V5),][,1:4]
length(share_snps$V1)
share_snps <- df2[(df2$V5 %in% df1$V5),][,1:4]
length(share_snps$V1)

df1_pt <- df1[!(df1$V5 %in% df2$V5),][,1:4]
length(df1_pt$V1)
df2_pt <- df2[!(df2$V5 %in% df1$V5),][,1:4]
length(df2_pt$V1)
write.table(share_snps, "./Del19_20_0p7x_share_snps.list", sep="\t", quote=F, row.names=F, col.names=F)
#write.table(df1_pt, "./Del19_private_site.list", sep="\t", quote=F, row.names=F, col.names=F)
#write.table(df2_pt, "./Del20_private_site.list", sep="\t", quote=F, row.names=F, col.names=F)
#write.table(share_snps, "./Del19_20_global_share_snps.list", sep="\t", quote=F, row.names=F, col.names=F)

find_shared <- function(f1_name,f2_name) {
  file1 = f1_name
  df1 <- read.delim(file1, header = FALSE, sep='\t')
  df1$V5 <- paste(df1$V1, df1$V2, sep="_")
  message("File ", f1_name, " SNPs: ", length(df1$V5))
  file2 = f2_name
  df2 <- read.delim(file2, header = FALSE, sep='\t')
  df2$V5 <- paste(df2$V1, df2$V2, sep="_")
  message("File ", f2_name, " SNPs: ", length(df2$V5))
  share_snps <- df1[(df1$V5 %in% df2$V5),][,1:4]
  message("Shared SNPs: ", length(share_snps$V1))
  df1_pt <- df1[!(df1$V5 %in% df2$V5),][,1:4]
  message("Private SNPs in ", f1_name, ": ", length(df1_pt$V1)) 
  #write.table(df1_pt, "./Del19_20_1p5x_private_snps.list", sep="\t", quote=F, row.names=F, col.names=F)
  df2_pt <- df2[!(df2$V5 %in% df1$V5),][,1:4]
  message("Private SNPs in ", f2_name, ": ", length(df2_pt$V1))  
}

find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "Del20_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_0p5x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_0p6x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_0p7x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_0p8x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_1x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_1p2x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")
find_shared("Del19_challenge_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", "down_sampling_1p5x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked.snplist.txt")


find_shared("Del19_20_global_share_snps.list", "Del19_20_1x_share_snps.list")
find_shared("Del19_20_1p5x_share_snps.list", "Del19_20_1p2x_share_snps.list")
