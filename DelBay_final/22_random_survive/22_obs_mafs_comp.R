setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/trash/obs_results_challenge_20")
df1 = read.delim("Challenge_20_all_minq20_minmq30_CV30_masked.mafs.test", header = TRUE, sep='\t')
head(df1)
df7 = read.delim("Sur_20_random_survive_minq20_minmq30_CV30_masked.mafs.test", header = TRUE, sep='\t')
df8 = read.delim("Ref_20_random_survive_minq20_minmq30_CV30_masked.mafs.test", header = TRUE, sep='\t')
     
df2 = read.delim("dfile_lrt_obs_test", header = FALSE, sep=' ')
head(df2)

df3 = read.delim("dfile_lrt_unkwn", header = FALSE, sep=' ')
head(df3)

df9 = read.delim("Challenge_20_random_survive_minq20_minmq30_CV30_masked_testassoc.mafs", header = TRUE, sep='\t')
df10 = read.delim("Challenge_20_random_survive_minq20_minmq30_CV30_masked.mafs", header = TRUE, sep='\t')

cor(df9$knownEM, df10$knownEM) 
plot(df9$knownEM, df10$knownEM)

cor(df8$knownEM, df2$V2) 
plot(df8$knownEM, df2$V2)

cor(df7$knownEM, df2$V3) 
plot(df7$knownEM, df2$V3)

cor(df1$knownEM, df2$V1) 
plot(df1$knownEM, df2$V1)

cor(df1$knownEM, df2$V2)          
plot(df1$knownEM, df2$V2)

cor(df2$V2, df3$V2)          
plot(df2$V2, df3$V2)

dim(df1[which(abs(df1$knownEM - df2$V2)>0.1),])


###########################
# Prepare global SNP list #
###########################
library(dplyr)
library(stringr)
setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/22_random_survive/Wild_18/")
df_mafs = read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", header = FALSE, sep='\t')
df_mafs = df_mafs[
  order( df_mafs[,1], df_mafs[,2] ),
]

write.table(df_mafs, file = "All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.sorted.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


df_testasso = read.delim("dfile_obsmaps", header = FALSE, sep='\t')

df_testasso %>% 
  mutate(V2 = str_replace(V2, "major=0", "A")) %>% 
  mutate(V2 = str_replace(V2, "major=1", "C")) %>% 
  mutate(V2 = str_replace(V2, "major=2", "G")) %>% 
  mutate(V2 = str_replace(V2, "major=3", "T")) %>% 
  mutate(V3 = str_replace(V3, "minor=0", "A")) %>% 
  mutate(V3 = str_replace(V3, "minor=1", "C")) %>% 
  mutate(V3 = str_replace(V3, "minor=2", "G")) %>% 
  mutate(V3 = str_replace(V3, "minor=3", "T"))  -> df_testasso

# check for major minor difference
# matched major - major and minor - minor
sum(df_mafs$V3 == df_testasso$V2 & df_mafs$V4 == df_testasso$V3)
# 2085137
sum(df_mafs$V3 == df_testasso$V2 & df_testasso$V3 == df_testasso$V2)
# 733
a<-df_testasso[df_mafs$V3 == df_testasso$V2 & df_mafs$V4 != df_testasso$V3, ]
sum(a$V2 == a$V3)
# 733
sum(a$V2 != a$V3)
# 385
sum(df_mafs$V3 != df_testasso$V2)
# 43904
# major = minor & minor = major
sum(df_mafs$V3 == df_testasso$V3 & df_mafs$V4 == df_testasso$V2)
# 43904

final_df = data.frame(df_mafs$V1, df_mafs$V2, df_testasso$V2, df_testasso$V3)
dim(final_df)
# swap the major and minor alleles if mismatches exist
final_df[df_testasso$V2 != df_mafs$V3, c("df_testasso.V2", "df_testasso.V3")] <- final_df[df_testasso$V2 != df_mafs$V3, c("df_testasso.V3", "df_testasso.V2")]

write.table(final_df, file = "Wild21_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs_testassoc_correction.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

