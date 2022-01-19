library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(related)
# to install export, one need to isntall X11(https://www.xquartz.org/), and then type devtools::install_github("tomwenseleers/export")
library(export)
library(qvalue)
library(stringr)
source("manhattan.R")

##############################
# check thinned pca results ##
##############################

setwd("~/Documents/HG/Domestication/04_pcadapt/")
library("pcadapt")
# pruned SNP set with 10K as the pruning window
path_to_file <- "./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("10K_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
graph2ppt(file="screeplot.pptx", width=5, height=4)
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
sample_file = read.delim("pop_509_sample_list.txt", header = TRUE, sep='\t')
poplist.names <- sample_file$Pop
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 10)
#########    DBW1      DBW2         DBX1      DBX2         DBX3      LIW1        LIW2        MEH2       MEW1      MEW2         NCW1      NCW2         NEH1      NEH2       UMFS       UNC1       UNC2
col <- c( "#33A65C", "#57A337",  "#F64971", "#FC719E", "#EB73B3",  "#30BCAD", "#21B087",  "#F89217", "#1BA3C6", "#2CB5C0",  "#A2B627", "#F8B620", "#F06719", "#E03426",  "#F8B620", "#CE69BE", "#A26DC2")

#col <- c("#26497a", "#904994", "#ea426d", "#ff7d00", "#00bf0d", "#82de86", "#86bde2")
#jpeg("Ind509_10K_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names, col=col)
graph2ppt(file="Ind509_10K_PC1-2.pptx", width=5, height=4)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names, col=col)
graph2ppt(file="Ind509_10K_PC2-3.pptx", width=5, height=4)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names, col=col)
graph2ppt(file="Ind509_10K_PC1-3.pptx", width=5, height=4)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names, col=col)
graph2ppt(file="Ind509_10K_PC1-4.pptx", width=5, height=4)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names, col=col)
graph2ppt(file="Ind509_10K_PC1-5.pptx", width=5, height=4)
dev.off()

x <- pcadapt(file, K = 2) # analysis with K = 2
# make PCA plot for thinned data
jpeg("Ind509_thinned_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names, col=col)
dev.off()
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
result_p <- x$pvalues
# adjust p-values for thinned dataset
qval <- qvalue(result_p)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers) # n = 276
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=paste0(CHR,'_',POS), Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Ind509_10K_thin_manhattan_by_chr.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 15),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Domestication n=510 K=2 outlier=276",)
dev.off()

##############################
###  start pcadapt naive #####
##############################

setwd("~/Documents/HG/Domestication/04_pcadapt/")

# pruned SNP set with 10K as the pruning window
path_to_file <- "./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("Ind509_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
sample_file = read.delim("pop_509_sample_list.txt", header = TRUE, sep='\t')
poplist.names <- sample_file$Pop
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("Ind509_naive_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names, col=col)
dev.off()
graph2ppt(file="Ind509_naive_PC1-2.pptx", width=5, height=3)

# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
jpeg("Ind509_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x , option = "manhattan")
dev.off()
# adjust p-values for naive dataset
result_p <- x$pvalues
qval <- qvalue(result_p)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=paste0(CHR,'_',POS), Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Ind509_naive_manhattan_by_chr.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Domestication naive n=510 K=2 outlier=807",)
dev.off()

##############################
###### start pcadapt BP ######
##############################

setwd("~/Documents/HG/Domestication/04_pcadapt/")
toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

#remove.packages("bigsnpr")
#install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
library(bigsnpr)
setwd("~/Documents/HG/Domestication/04_pcadapt/")

# all
bedfile = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# obtain the "bed" snp index during pruning and manually store them into a txt file
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8)
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = "all_pruned_SNP_list.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2]))
# load the na_idx
#na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
#Gm = toMatrix(G)
#Gm = Gm[,-na_idx]
#G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
#test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
hist(result_p, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
# Choosing a cutoff for outlier detection
qval <- qvalue(result_p)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)

# FDR correction
#padj <- p.adjust(result_p,method="BH")
#alpha <- 0.05
#outliers <- which(padj < alpha)
#length(outliers)

print(CHR[outliers])
daf = data.frame(CHR=CHR, POS=POS, SNP=paste0(CHR,"_",POS), Ps=result_p, qvalue=qval)
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
outlier_daf = data.frame(CHR=CHR[outliers], POS=POS[outliers], POS=POS[outliers], snp_id=outlier_SNP)
write.table(outlier_daf, file = "PCAdapt_outliers_q05_n_428.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_Ind509_best_practice_outlier_PC1-2_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="qvalue", highlight1=outlier_SNP, logp=T, cex.axis = 0.8, ylim = c(0, 6),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(q-value)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()



