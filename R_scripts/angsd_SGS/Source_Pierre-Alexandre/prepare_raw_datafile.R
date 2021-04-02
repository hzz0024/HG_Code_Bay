# setting the working directory
setwd ("C:/Users/Pierre-Alexandre/Dropbox/RAD_Dorade_Anchois_Sole/Post_Stacks/DATA_batch_2&3/GDS")

## Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)

genofile <- openfn.gds("gds_batch_2_couv_hwe_ok_snp.gds")

## make a PED file for PLINK
snpset <- snpgdsSelectSNP(genofile, verbose=TRUE, autosome.only=FALSE)
# save the list of SNPs to a file
write.table(snpset,"snpset_34679.txt",row.names = FALSE,col.names = FALSE)
# Convert to PLINK format to produce a PED file and a MAP file
snpgdsGDS2PED(genofile, ped.fn="batch_2_couv_hwe_ok_snp", snp.id=snpset)
