library("pcadapt")
library(export)
library(stringr)
source("manhattan.R")
# all
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("all_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
jpeg("PCloading_full_PC1-6.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs loadings PC", 6))
dev.off()


# 1K
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

jpeg("PCloading_1K_PC1-6.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned loadings PC", 6))
dev.off()

poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
#jpeg("1K_thin_pca.jpg", width = 16, height = 9, units = 'in', res = 300)
#par(mfrow=c(3,2))
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="1K_PC1-2.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="1K_PC1-3.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="1K_PC1-4.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="1K_PC1-5.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 6, pop = poplist.names)
graph2ppt(file="1K_PC1-6.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 4, j = 5, pop = poplist.names)
graph2ppt(file="1K_PC4-5.pptx", width=4.5, height=2.5)
#dev.off()


# all no inversions
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("all_noinvers_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

jpeg("PCloading_full_noinvers_PC1-6.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 6))
dev.off()

# all no inversions no chr56 inversions
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_nochr56invers.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_nochr56invers.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("all_noinvers_nochr56invers_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

jpeg("PCloading_full_noinvers_nochr56invers_PC1-6.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.002, 0.002),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("Full SNPs no invers loadings PC", 6))
dev.off()



# 1K no invers
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.thinned.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.thinned.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 5)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("1K_noinvers_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

jpeg("test.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers loadings PC", 6))
dev.off()

# 1K no invers no chr56 invers
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.nochr56invers.thinned.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.nochr56invers.thinned.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("1K_noinvers_nochr56invers_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# check loadings
#par(mfrow = c(2, 2))
#for (i in 1:4)
#  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

jpeg("PCloading_1K_noinvers_nochr56invers_PC1-6.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 1])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 1))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 2])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 2))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 3])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 3))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 4])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 4))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 5])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 5))
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ls=x$loadings[, 6])
manhattan(daf, chr="CHR",bp="POS", p="Ls", logp=FALSE, cex.axis = 0.8, ylim = c(-0.005, 0.005),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="loadings", cex.lab=1.4, main =paste0("1K thinned no invers & chunked chr5-6 loadings PC", 6))
dev.off()

poplist.names <- c(rep("DelBay_CS", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_Sel_NEH", 6), rep("Louisiana_Sel", 6), rep("Louisiana_SL", 6))
print(poplist.names)
# Computing the test statistic based on PCA
#jpeg("1K_thin_pca.jpg", width = 16, height = 9, units = 'in', res = 300)
#par(mfrow=c(3,2))
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC1-2.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC1-3.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC1-4.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC1-5.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 1, j = 6, pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC1-6.pptx", width=4.5, height=2.5)
plot(x, option = "scores", i = 4, j = 5, pop = poplist.names)
graph2ppt(file="1K_noinverschr56inves_PC4-5.pptx", width=4.5, height=2.5)
#dev.off()
