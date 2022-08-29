#The minimal number of SNPs in a ROH was
#determined by the formula proposed by Lencz et al.
#and adapted by Purfield et al.

a = 0.05  # percentage of false positive ROH
ns = 141960  # number of genotyped SNPs per individual
ni = 509 # the number of genotyped individuals
het = 0.2413 # het the mean heterozygosity across all SNPs

L = (log(a/ns*ni))/log(1-het)
L

# plink.hom
# 