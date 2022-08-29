
###############################################
# Data preparation
# 1. convert the vcf or plink format files to GENEPOP format using PGDspider
# 2. change the last sample ID in each population to its population name, also useful for the DAPC figure legend plotting

##### for example, change 
##### Pop
##### WEE.10,	0102	0404	0202	0101	0101
##### WEE.11,	0101	0404	0202	0101	0101
##### .....
##### WEE.23,	0101	0404	0202	0202	0101

##### into (only change WEE.23 to HAT here)

##### Pop
##### WEE.10,	0102	0404	0202	0101	0101
##### WEE.11,	0101	0404	0202	0101	0101
##### .....
##### HAT,	0101	0404	0202	0202	0101
##### here "HAT" is the population name for this group

# note: this example genepop file is converted from SNP genotyping data from MassARRAY system, PGDspider can do the formatting but need to specify a population definition file
###############################################

library(adegenet)
library(Demerelate)
library(dartR)
library(adegenet)
library(vcfR)
library(export)
adegenetServer(what = "DAPC")

setwd("~/Dropbox/Mac/Documents/HG/Domestication/10_DAPC_outlier")

vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --bed pop_n_477_pcadapt_outflank.shared.10k.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv --recode --recode-INFO-all --out n_509_shared_outliers", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf
# --recode-INFO-all
# --out n_509_shared_outliers
# --recode
# --bed pop_n_477_pcadapt_outflank.shared.10k.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv
# 
# After filtering, kept 509 out of 509 Individuals
# Outputting VCF file...
# Read 35 BED file entries.
# After filtering, kept 410 out of a possible 141960 Sites
# Run Time = 2.00 seconds

# load the population information
pop_info <- read.table("pop_509_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
# load vcf file
vcf_file = "n_509_shared_outliers.recode.vcf"
#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
x <- Mydata1

sum(is.na(Mydata1$tab))

# estimate the number of PCs for data analysis
dapc1 <- dapc(x, x@pop)
temp <- optim.a.score(dapc1)
library(export)
graph2ppt(file="Domestication_alpha_score",width=12,height=12)

# temp$best = 26
x@pop
x@all.names
grp <- find.clusters(Mydata1, max.n.clust=20)
names(grp)
head(grp$grp, 10)
table(pop(x), grp$grp)
table.value(table(pop(x), grp$grp), col.lab=levels(pop(x)),
            row.lab=paste("ori", 1:17))
# formal running
dapc1 <- dapc(x, x@pop, n.pca = 26, n.da=10)
scatter(dapc1)
# adegenetServer(what = "DAPC")
# #reading the genepop file as input
# x <- read.genepop("WAL62_SNP1978.gen",quiet = TRUE)
# x <- read.genepop("WAL605_SNP68.gen",quiet = TRUE)
# x <- read.genepop("WAL547_SNP68.gen",quiet = TRUE)

#dapc1 <- dapc(x, var.contrib=FALSE, scale=FALSE, n.pca=50, n.da=1)

#plot
# cstar=0 will remove the connection lines
# add cell = 0  will remove the ellipses
# cex is the size of the dot
#SMB plot

#   MEW1       MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
col1 <- c("#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA", 
               #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
               "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")

col_transparent <- adjustcolor(col1,alpha.f = 0.1)

scatter(dapc1,  xlab="Discriminatn Function 1", ylab="Discriminatn Function", cex = 1.5, pch=c(rep(18, 9), rep(20,8)), lwd=2, lty=2, solid = .8, cstar=0, scree.da=FALSE,  col = col1, clab = 1 ) #label.inds = list(air = 2, pch = NA))
legend("right", legend = levels(pop(x)),cex = 1,col = col1, bty = 'n', pch=c(rep(18, 9), rep(20,8)), xpd = TRUE, inset=c(0.2,0.55))
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=15, type="h", lwd=3)
}
add.scatter(myInset(), posi="topleft",
            inset=c(0.05,0.55), ratio=.14,
            bg=transp("white"))

library(export)
graph2ppt(file="Domestication_outlier_DAPC",width=8,height=9)
####################################################


