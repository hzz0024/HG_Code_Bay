# this R script is used to create the prunning snp list for each pairwise population

setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt/")
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_DEBY.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_DEBY.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "CS_DEBY_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_NEH.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_NEH.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "CS_NEH_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.SL_OBOYS2.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.SL_OBOYS2.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "SL_OBOYS2_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "all_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

library(bigsnpr)
setwd("/Volumes/cornell/CVreseq_random_forest")
bedfile = "CS_HC.chr2.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("CS_HC.chr2.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
#newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, size = 10)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "CS_HC.chr2.clumping.list", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)


bedfile = "HCVA_CLP.chr2.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("HCVA_CLP.chr2.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
#newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, size = 10)
which_pruned =  attr(newpc, which="subset")
SNP_list <- SNPs[which_pruned]
write.table(SNP_list, file = "HCVA_CLP.chr2.clumping.list", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
