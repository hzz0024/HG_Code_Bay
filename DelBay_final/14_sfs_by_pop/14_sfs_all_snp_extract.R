########## format for global SNP list ###########

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/21_ngsparalog")

# below is the code to extract non-paralog SNPs (global snp with p-value and maf filtering)
All <- read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.txt", header = FALSE, sep='\t')
head(All)
All$snp <- paste0(All$V1, "_", All$V2)
dim(All)

non_paralogs <- read.delim("snps_paralogs_fdr05.txt", header = FALSE, sep=' ')
non_paralogs$snp <- paste0(non_paralogs$V1, "_", non_paralogs$V2)
dim(non_paralogs)
non_paralogs <- All[-which(All$snp %in% non_paralogs$snp),]
dim(non_paralogs)
non_paralogs <- non_paralogs[with(non_paralogs, order(V1, V2)), ]
write.table(non_paralogs[,1:4],"All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_noparalogs.snplist.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

# This is not the END, we need to adjust the major and minor alleles based on testassoc analysis. See step 22!!!!!!!!!!!!!!!

########## format for sfs SNP list ###########

setwd("~/Dropbox/Mac/Documents/HG/DelBay_final/14_SFS")
# below is the code to extract non-paralog SNPs from sfs sites
All <- read.delim("All_sfs_all_sites.snplist.sorted.txt", header = FALSE, sep='\t')
head(All)

All$snp <- paste0(All$V1, "_", All$V2)

paralogs <- read.delim("snps_paralogs_fdr05.txt", header = FALSE, sep=' ')
paralogs$snp <- paste0(paralogs$V1, "_", paralogs$V2)
non_paralogs <- All[-which(All$snp %in% paralogs$snp),]
dim(non_paralogs)
#[1] 30964046        5
non_paralogs <- non_paralogs[with(non_paralogs, order(V1, V2)), ]
write.table(non_paralogs[,1:4],"All_sfs_all_sites_noparalogs.snplist.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


# outliers from SGS and permutation
SGS_wild_18 <- read.delim("18_SGS_HC_NB_ps_outlier.list", header = TRUE, sep='\t')$id
SGS_wild_19 <- read.delim("19_SGS_HC_NB_ps_outlier.list", header = TRUE, sep='\t')$id
SGS_wild_21 <- read.delim("21_SGS_HC_NB_ps_outlier.list", header = TRUE, sep='\t')$id
SGS_challenge_19 <- read.delim("19_SGS_Sur_Ref_ps_outlier.list", header = TRUE, sep='\t')$id
SGS_challenge_20 <- read.delim("20_SGS_Sur_Ref_ps_outlier.list", header = TRUE, sep='\t')$id
perm_wild_18 <- read.delim("18_HC_NB_permutation_outlier.txt", header = TRUE, sep='\t')$SNP
perm_wild_19 <- read.delim("19_HC_NB_permutation_outlier.txt", header = TRUE, sep='\t')$SNP
perm_wild_21 <- read.delim("21_HC_NB_permutation_outlier.txt", header = TRUE, sep='\t')$SNP
perm_challenge_19 <- read.delim("19_Sur_Ref_permutation_outlier.txt", header = TRUE, sep='\t')$SNP
perm_challenge_20 <- read.delim("20_Sur_Ref_permutation_outlier.txt", header = TRUE, sep='\t')$SNP

outlier_list <- c(SGS_wild_18, SGS_wild_19, SGS_wild_21, SGS_challenge_19, SGS_challenge_20, perm_wild_18, perm_wild_19, perm_wild_21,  perm_challenge_19, perm_challenge_20)

non_paralogs_outliers <- non_paralogs[-which(non_paralogs$snp %in% outlier_list),]
dim(non_paralogs_outliers)
#[1] 30904801        5
non_paralogs_outliers <- non_paralogs_outliers[with(non_paralogs_outliers, order(V1, V2)), ]

# random sample 1M SNPs
random_1M <- sample_n(non_paralogs_outliers, 1000000)
random_1M <- random_1M[with(random_1M, order(V1, V2)), ]
write.table(random_1M[,1:4],"All_sfs_all_sites_noparalogs_random_1M.snplist.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

########## plot for LR ###########
lr <- read.table("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.snplist.output")
lr$pval <- 0.5*pchisq(lr$V5,df=1,lower.tail=FALSE) # append column of p-values
lr$pval.adj <- p.adjust(lr$pval, method="BH") # p-values adjusted for number of tested sites
# The 7th column of the lr data.frame is the adjusted p-value for rejecting the null hypothesis that reads
# covering the site derive from a single locus. Of course you can use any p-value adjustment of your
# choosing, e.g. "fdr".

lr %>% 
  mutate(V1 = str_replace(V1, "NC_035780.1", "1")) %>% 
  mutate(V1 = str_replace(V1, "NC_035781.1", "2")) %>% 
  mutate(V1 = str_replace(V1, "NC_035782.1", "3")) %>% 
  mutate(V1 = str_replace(V1, "NC_035783.1", "4")) %>% 
  mutate(V1 = str_replace(V1, "NC_035784.1", "5")) %>% 
  mutate(V1 = str_replace(V1, "NC_035785.1", "6")) %>% 
  mutate(V1 = str_replace(V1, "NC_035786.1", "7")) %>% 
  mutate(V1 = str_replace(V1, "NC_035787.1", "8")) %>%
  mutate(V1 = str_replace(V1, "NC_035788.1", "9")) %>% 
  mutate(V1 = str_replace(V1, "NC_035789.1", "10"))  -> lr
lr$V1 <- as.numeric(lr$V1) 
head(lr)
lr$SNP <- paste0(lr$V1, "_", lr$V2)
lr$pval.adj[lr$pval.adj == 0] = 7.539667e-309
source("manhattan.R")
manhattan(lr,chr="V1",bp="V2",p="pval.adj",logp=TRUE, ylim = c(0, max(-log10(lr$pval.adj))),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=T,
          ylab="Mismapping LR") 
