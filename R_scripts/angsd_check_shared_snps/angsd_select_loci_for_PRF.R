setwd("~/Documents/Ryan_workplace/DelBay_adult/batch_effect/downsampling")
#target <- read.delim("PRF_toy_Del20_FDR_outlier_rd1.list", header = TRUE, sep='\t')
#target$id = paste0(target$chr,'_',target$pos)
#head(target)

df1 <- read.delim("Del20_final_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.mafs", header = TRUE, sep='\t')
df1$id = paste0(df1$chromo,'_',df1$position)
length(df1$chromo)

df2 <- read.delim("Del19_20_1x_share_snps.list", header = FALSE, sep='\t')
df2$id = paste0(df2$V1,'_',df2$V2)
length(df2$id)

df3 <- df1[(df1$id %in% df2$id),]
df3 <- df3[with(df3, order(knownEM)),]
length(df3$chromo)

rd1 <- df3[1:182,]
write.table(rd1, "random182_maf05.txt", row.names = FALSE, sep="\t", quote = FALSE)

# start to extract maf information for target SNPs
target_snps <- df2[(df2$id %in% target$id),]
target_snps
