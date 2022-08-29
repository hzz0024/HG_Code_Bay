library("readr")
library("ggplot2")
library("dplyr")
require(gridExtra)
devtools::install_github('Mikata-Project/ggthemr')
library(ggthemr)
# https://speciationgenomics.github.io/filtering_vcfs/ for detailed steps
setwd("~/Documents/HG/Domestication/00_vcf")
################################
###### SNP missingness #########
################################
library(ggthemr)
ggthemr("fresh")
#ggthemr_reset()
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-indv --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe", sep=""))
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-site --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe", sep=""))
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --freq2 --max-alleles 2 --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe", sep=""))

var_miss <- read_delim("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_miss

a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
aa <- a + ylab("Density") + xlab("proportion of missing data per SNP")
summary(var_miss$fmiss)
################################
###### ind missingness #########
################################
ind_miss  <- read_delim("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
b <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
bb<-b + ylab("Sample count") + xlab("proportion of missing data per individual") 
##############################################
###### allele frequency distribution #########
##############################################
var_freq <- read_delim("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
c <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
cc<-c + ylab("Density") + xlab("Minor allele frequency (MAF)") 

jpeg("vcf_summary.jpg", width = 6, height = 8, units = 'in', res = 300)

grid.arrange(aa, bb, cc, nrow=3)
dev.off()



