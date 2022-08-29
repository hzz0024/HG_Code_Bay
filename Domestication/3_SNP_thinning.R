# this R script is used to create the prunning snp list for each pairwise population
library(bigsnpr)

# set up the directory
setwd("~/Documents/HG/Domestication/03_LD_clumping/")

# define toMatrix function to create a matrix
toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

sub_name = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral"

f_bk = paste0(sub_name, ".bk")
if (file.exists(f_bk)) {
  #Delete file if it exists
  file.remove(f_bk)
}

#  data preparation
snp_readBed(paste0(sub_name, ".bed"))
# this will create a .rds file
obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# genotype imputation
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
#big_counts(G, ind.col = 1:dim(G)[1]) # check if NAs are 0
# LD clumping using r2 = 0.2 and window size of 10K
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) # size of LD clumping window is 10K
# extract SNPs after clumpping
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = paste0(sub_name, "_LDclump_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(paste0("SNPs after clumpping is: ", length(keep_snp_ids), " out of ", dim(obj.bigSNP$map)[1]))
# generate LD clumped vcf file
system(paste(vcftools," --vcf ",sub_name,".recode.vcf", " --snps ", sub_name, "_LDclump_SNP.txt", " --recode --recode-INFO-all --out ", sub_name, "_pruned", sep=""))
#system(paste(plink, " --vcf ",sub_name,"_pruned.recode.vcf", " --allow-extra-chr --make-bed --out ", sub_name,"_pruned", sep=""))



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
