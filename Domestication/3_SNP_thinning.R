# this R script is used to create the prunning snp list for each pairwise population
library(bigsnpr)

# set up the directory
setwd("~/Documents/HG/Domestication/03_SNP_thinning")

bedfile = "genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("genetyped_data_n_514_maf05_maxmiss095_popmiss095_hwe.rds")
# impute the genotypes based on snp_fastImputeSimple and method (see https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html for details)
# Either "random" (sampling according to allele frequencies), "mean0" (rounded mean), "mean2" (rounded mean to 2 decimal places), "mode" (most frequent call).
G <- snp_fastImputeSimple(obj.bigSNP$genotypes, method = c("mean0"), ncores = 8)
#G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, size = 10)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "all_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
