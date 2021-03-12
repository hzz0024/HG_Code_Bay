# this script is used to find common shared snps between two snp list files
file1 = "Del19_final_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers.snplist.txt"
df1 <- read.delim(file1, header = FALSE, sep='\t')
df1$V5 <- paste(df1$V1, df1$V2, sep="_")
file2 = "Del20_final_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers.snplist.txt"
df2 <- read.delim(file2, header = FALSE, sep='\t')
df2$V5 <- paste(df2$V1, df2$V2, sep="_")

share_snps <- df1[(df1$V5 %in% df2$V5),][,1:4]
length(share_snps$V1)
share_snps <- df2[(df2$V5 %in% df1$V5),][,1:4]
length(share_snps$V1)

df1_pt <- df1[!(df1$V5 %in% df2$V5),][,1:4]
length(df1_pt$V1)
df2_pt <- df2[!(df2$V5 %in% df1$V5),][,1:4]
length(df2_pt$V1)
write.table(df1_pt, "./Del19_private_site.list", sep="\t", quote=F, row.names=F, col.names=F)
write.table(df2_pt, "./Del20_private_site.list", sep="\t", quote=F, row.names=F, col.names=F)
write.table(share_snps, "./Del19_20_global_share_snps.list", sep="\t", quote=F, row.names=F, col.names=F)

