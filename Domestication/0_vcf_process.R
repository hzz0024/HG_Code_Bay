library("readr")
library("ggplot2")
library("dplyr")
require(gridExtra)
# https://speciationgenomics.github.io/filtering_vcfs/ for detailed steps
setwd("~/Documents/HG/Domestication/00_vcf")
################################
###### SNP missingness #########
################################
# vcftools --vcf genetyped_data_n_514_maf05_maxmiss07_nochr156invers.recode.vcf --missing-indv
# vcftools --vcf genetyped_data_n_514_maf05_maxmiss07_nochr156invers.recode.vcf --missing-site
# vcftools --vcf genetyped_data_n_514_maf05_maxmiss07_nochr156invers.recode.vcf --freq2 --max-alleles 2
var_miss <- read_delim("./out.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_miss

a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
aa <- a + theme_light() + ylab("Density") + xlab("proportion of missing data per SNP")
summary(var_miss$fmiss)
################################
###### ind missingness #########
################################
ind_miss  <- read_delim("./out.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
b <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
bb<-b + theme_light() + ylab("Count") + xlab("proportion of missing data per individual") 
##############################################
###### allele frequency distribution #########
##############################################
var_freq <- read_delim("./out.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
c <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
cc<-c + theme_light() + xlab("Minor allele frequency (MAF)") 

jpeg("vcf_summary.jpg", width = 6, height = 8, units = 'in', res = 300)
grid.arrange(aa, bb, cc, nrow=3)
dev.off()



