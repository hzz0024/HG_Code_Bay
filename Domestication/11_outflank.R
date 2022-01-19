devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(stringr)
library(related)
library(bigstatsr)
library(vcfR)
library(RColorBrewer)
library(fields)
library(pcadapt)
library(devtools)
source("manhattan.R")
# http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html
# set up the working directory
setwd("~/Dropbox/Mac/Documents/HG/Domestication/11_outflank")

# define the Matrix function
toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

# load the data
bedfile = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- paste0(CHR, "_", POS)

# imputate the genotype "mean0" (rounded mean) https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html 
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8)
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
length(which_pruned)
keep_snp_ids = SNPs[which_pruned]
# start to calculate fst
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
# load the population information
pop <- read.delim("pop_n_509.txt",header = F)$V1
my_fst <- MakeDiploidFSTMat(toMatrix(G_coded), locusNames = paste0(CHR,'_',POS), popNames = pop)
# Data checks: Heterozygosity vs. FST
plot(my_fst$He, my_fst$FST)
plot(as.numeric(my_fst$LocusName), my_fst$He, ylim=c(0,0.5), pch=19, col=rgb(0,0,0,0.1))
# check the chi-square distributed when all loci are included 
hist(my_fst$FSTNoCorr, breaks=seq(0,1, by=0.01))
# Removing low Heterozygosity variants results in a more chi-square looking FST distribution
hist(my_fst$FSTNoCorr[my_fst$He>0.1], breaks=seq(0,1, by=0.01))
# Data checks: FST vs. FSTNoCorr
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)

### All data
k <- 17 ## Number of pops 
out_ini <- OutFLANK(my_fst, NumberOfSamples=k) 
# Plot results to compare chi-squared distribution vs. actual FST distribution
OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Poor fit, particularly on right tail
OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.05, titletext = NULL)
# Histogram of P-values 
hist(out_ini$results$pvaluesRightTail,breaks = 20)

### With LD pruning
#### Evaluating OutFLANK with pruned data ####
plot(my_fst$He[which_pruned], my_fst$FST[which_pruned])
Fstdf2 <- my_fst[which_pruned,] 
dim(Fstdf2)
Fstdf3 <- Fstdf2[Fstdf2$He>0.1,]
plot(Fstdf3$He, Fstdf3$FST)
### Trimming without He trimming
out_trim1 <- OutFLANK(Fstdf2, NumberOfSamples=k)
OutFLANKResultsPlotter(out_trim1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.15, titletext = NULL)
# Trimming loci with low He
out_trim <- OutFLANK(Fstdf3, NumberOfSamples=k)
head(out_trim$results)
plot(out_trim$results$He, out_trim$results$FST)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.15, titletext = NULL)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.05, titletext = NULL)


# best practice
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)
# retain the SNPs with He > Hmin (0.1)
Hmin = 0.1
keepers = which((P1$FSTNoCorr > 0) & (P1$He >= Hmin))
P1 = P1[keepers,]
# count how many outlier candidates with q-value < 0.05
length(P1$OutlierFlag[P1$OutlierFlag==TRUE])

# notice how the output is ordered differently
my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")

# plot the p-values
hist(P1$pvaluesRightTail)
# plot the q-values
hist(P1$qvalues)
# plot the qvalues using Mahattan plot
jpeg("Mahattan_Ind509_best_practice_outlier.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
CHR = sapply(strsplit(P1$LocusName, "_"), "[[", 1)
POS = sapply(strsplit(P1$LocusName, "_"), "[[", 2)
P1$CHR = as.numeric(CHR)
P1$POS = as.numeric(POS)
P1$SNP = P1$LocusName
manhattan(P1, chr="CHR",bp="POS", p="pvalues", highlight1=P1$LocusName[P1$OutlierFlag ==TRUE], logp=T, cex.axis = 0.8, ylim = c(0, 12),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(q-values)", cex.lab=1.4, main = "Outflank Best Practice")
dev.off()

Outlier = P1$LocusName[P1$OutlierFlag ==TRUE]
Outlier_chr = P1$CHR[P1$OutlierFlag ==TRUE]
Outlier_pos = P1$POS[P1$OutlierFlag ==TRUE]
df = data.frame(Outlier_chr, Outlier_pos, Outlier_pos, Outlier)

write.table(df, file = "Outflank_outliers_q05_n_1366.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
s###########################################################
# Bonus store pruned data with different thinning windows #
###########################################################

bedfile = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed"
# this will create a .rds file
snp_readBed(bedfile)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# Note that most of the algorithms of this package don’t handle missing values. 
# I used snp_fastImputeSimple() to impute missing values of genotyped variants.
#G2 <- snp_fastImputeSimple(G, method = c("random"))
#G3 <- snp_fastImputeSimple(G, method = c("mode"))
SIZE <- c(20, 50, 100, 200, 500)
#SIZE <- c(50)
ind_keeps = list()
# perform snp_clumping on masked dataset at different window size (in KB)
for(i in seq(length(SIZE))){
  size_ = SIZE[i]
  ind.keep_ <- snp_clumping(
    G,
    infos.chr = CHR,
    infos.pos = POS,
    ind.row = rows_along(G),
    S = NULL,
    thr.r2 = 0.2,
    size = size_,
    exclude = NULL,
    ncores = 1
  )
  ind_keeps[[i]] = ind.keep_
}
# 1-5 the index list
index <- ind_keeps[[2]]
# output the index list
write.table(obj.bigSNP$map$marker.ID[index], file = "SNP_thinned.txt", sep = "\t",
            row.names = FALSE, quote = F, col.names=FALSE)

          