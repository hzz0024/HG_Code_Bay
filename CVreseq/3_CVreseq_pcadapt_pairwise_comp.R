#library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(related)
library(export)
library(stringr)
source("manhattan.R")
library(export)

####################### check thinned pca results #####################
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt/thinned_vcf")

# 500K
library("pcadapt")
path_to_file <- "./500K.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./500K.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("500K_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 5)
jpeg("500K_thin_pca_k3.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="500K_PC1-2.pptx", width=5, height=3)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
graph2ppt(file="500K_PC2-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="500K_PC1-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="500K_PC1-4.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="500K_PC1-5.pptx", width=5, height=3)
dev.off()

# 5K
library("pcadapt")
path_to_file <- "./5K.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./5K.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("5K_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 5)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="5K_PC1-2.pptx", width=5, height=3)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
graph2ppt(file="5K_PC2-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="5K_PC1-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="5K_PC1-4.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="5K_PC1-5.pptx", width=5, height=3)

# 1K
library("pcadapt")
path_to_file <- "./1K.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./1K.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("1K_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
jpeg("PCloading_1K_PC2.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main = "1K Thinned PC2 loadings")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 5)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="1K_PC1-2.pptx", width=5, height=3)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
graph2ppt(file="1K_PC2-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="1K_PC1-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="1K_PC1-4.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="1K_PC1-5.pptx", width=5, height=3)

# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_DEBY_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()


# CS_DEBY
library("pcadapt")
path_to_file <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.CS_DEBY.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.CS_DEBY.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("CS_DEBY_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("CS_DEBY_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_DEBY_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("CS_DEBY_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs Ches_Sel_DEBY K=2",)
dev.off()

# CS_NEH
library("pcadapt")
path_to_file <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.CS_NEH.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.CS_NEH.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("CS_NEH_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("DelBay_Sel_NEH", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("CS_NEH_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_NEH_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("CS_NEH_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs DelBay_Sel_NEH K=2",)
dev.off()


# SL_OBOYS2
library("pcadapt")
path_to_file <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.SL_OBOYS2.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./Thinned.SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.nomissing.FIL.format.SL_OBOYS2.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("SL_OBOYS2_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_Sel", 6),rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("SL_OBOYS2_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("SL_OBOYS2_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("SL_OBOYS2_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Louisiana_Sel vs Louisiana_SL K=2",)
dev.off()

# all five populations
library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.thin.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.thin.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("all_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
graph2ppt(file="all_thin_screeplot.pptx", width=7, height=5)
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 3)
jpeg("all_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="all_thin_pca_k2.pptx", width=7, height=5)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("SL_OBOYS2_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("all_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4)
dev.off()

############################ start pcadapt naive #########################
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt")
# CS_DEBY
library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_DEBY.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_DEBY.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("CS_DEBY_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("CS_DEBY_naive_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_DEBY_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("CS_DEBY_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 160),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs Ches_Sel_DEBY K=2",)
dev.off()

# CS_NEH
library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_NEH.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.CS_NEH.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("CS_NEH_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("DelBay_Sel_NEH", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("CS_NEH_naive_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_NEH_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("CS_NEH_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs DelBay_Sel_NEH K=2",)
dev.off()


# SL_OBOYS2
library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.SL_OBOYS2.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.maf05.nomissing.SL_OBOYS2.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("SL_OBOYS2_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_Sel", 6),rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("SL_OBOYS2_naive_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("SL_OBOYS2_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("SL_OBOYS2_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 50),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Louisiana_Sel vs Louisiana_SL K=2",)
dev.off()


# all
library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("all_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 3)
jpeg("all_pca_naive_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
#jpeg("CS_DEBY_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
#plot(x , option = "manhattan")
#dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("all_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 200),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4)
dev.off()

############################ start pcadapt BP #########################
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt")

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

remove.packages("bigsnpr")
install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
library(bigsnpr)
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt/")

# CS_DEBY
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
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("CS_DEBY_filter_idx.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
write.table(outlier_SNP, file = "CS_DEBY_outlier_SNP.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_CS_DEBY_BP.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs Ches_Sel_DEBY Best Practice")
dev.off()

# CS_NEH
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
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("CS_NEH_filter_idx.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
write.table(outlier_SNP, file = "CS_NEH_outlier_SNP.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_CS_NEH_BP.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "DelBay_CS vs DelBay_Sel_NEH Best Practice")
dev.off()

# SL_OBOYS2
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
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("SL_OBOYS2_filter_idx.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
write.table(outlier_SNP, file = "SL_OBOYS2_outlier_SNP.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_SL_OBOYS2_BP.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Louisiana_Sel vs Louisiana_SL Best Practice")
dev.off()


# all
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 1, min.mac = 1) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = "keep_snp_ids_1K.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
outlier_daf = data.frame(CHR=CHR[-na_idx][outliers], POS=POS[-na_idx][outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "all_BP_outlier_PC1-5_5K.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_all_BP_outlier_PC1-5_5K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()


############ all without oboys2_5 outlier ################
# all
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noOBOY5.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noOBOY5.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 1, min.mac = 1) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = "keep_snp_ids_noOBOYS5_1K.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_430.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
outlier_daf = data.frame(CHR=CHR[-na_idx][outliers], POS=POS[-na_idx][outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "all_BP_outlier_PC1-5_5K_noOBOYS5.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_all_BP_outlier_PC1-5_5K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()


# full no inversions
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

bedfile_noinvers = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.bed"
snp_readBed(bedfile_noinvers)
# this will create a .rds file
obj.bigSNP_noinvers <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.rds")
G_noinvers <- obj.bigSNP_noinvers$genotypes
SNPs_noinvers <- obj.bigSNP_noinvers$map$marker.ID
CHR_noinvers <- obj.bigSNP_noinvers$map$chromosome
POS_noinvers <- obj.bigSNP_noinvers$map$physical.pos

# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G_noinvers, infos.chr = CHR_noinvers, infos.pos = POS_noinvers, thr.r2 = 0.2, size = 1, min.mac = 1, k = 2) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs_noinvers[which_pruned]
write.table(keep_snp_ids, file = "keep_snp_ids_1K_noinvers.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
outlier_daf = data.frame(CHR=CHR[-na_idx][outliers], POS=POS[-na_idx][outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "all_BP_outlier_PC1-2_1k_noinvers.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_all_BP_outlier_PC1-2_1K_noinvers.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

# Full no inversions no chr56 inversions
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

bedfile_noinvers = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_nochr56invers.bed"
snp_readBed(bedfile_noinvers)
# this will create a .rds file
obj.bigSNP_noinvers <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_nochr56invers.rds")
G_noinvers <- obj.bigSNP_noinvers$genotypes
SNPs_noinvers <- obj.bigSNP_noinvers$map$marker.ID
CHR_noinvers <- obj.bigSNP_noinvers$map$chromosome
POS_noinvers <- obj.bigSNP_noinvers$map$physical.pos

# obtain the "bed" snp index during pruning and manually store them into a txt file
newpc <- snp_autoSVD(G_noinvers, infos.chr = CHR_noinvers, infos.pos = POS_noinvers, thr.r2 = 0.2, size = 500, min.mac = 1) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs_noinvers[which_pruned]
write.table(keep_snp_ids, file = "keep_snp_ids_1K_noinvers_nochr56invers.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
# load the na_idx
na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[-na_idx][outliers])
daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])
outlier_daf = data.frame(CHR=CHR[-na_idx][outliers], POS=POS[-na_idx][outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "all_BP_outlier_PC1-5_500k_noinvers_nochr56invers.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_all_BP_outlier_PC1-5_500K_noinvers_nochr56invers.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

