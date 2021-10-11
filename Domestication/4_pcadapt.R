library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(related)
# to install export, one need to isntall X11(https://www.xquartz.org/), and then type devtools::install_github("tomwenseleers/export")
library(export)
library(stringr)
source("manhattan.R")
####################### check thinned pca results #####################
setwd("~/Documents/HG/Domestication/04_pcadapt/maf005/")
library("pcadapt")
# pruned SNP set with 10K as the pruning window
path_to_file <- "./genetyped_data_n_514_maf05_maxmiss07_pruned.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./genetyped_data_n_514_maf05_maxmiss07_pruned.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("500K_thin_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
sample_file = read.delim("pop_514_sample_list.txt", header = TRUE, sep='\t')
poplist.names <- sample_file$Pop
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("Ind514_10K_thin_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file="Ind514_10K_PC1-2.pptx", width=5, height=3)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
graph2ppt(file="Ind514_10K_PC2-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 3, pop = poplist.names)
graph2ppt(file="Ind514_10K_PC1-3.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)
graph2ppt(file="Ind514_10K_PC1-4.pptx", width=5, height=3)
plot(x, option = "scores", i = 1, j = 5, pop = poplist.names)
graph2ppt(file="Ind514_10K_PC1-5.pptx", width=5, height=3)
dev.off()

# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
jpeg("Ind514_10K_thin_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x , option = "manhattan")
dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=paste0(CHR,'_',POS), Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Ind514_10K_thin_manhattan_by_chr.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 15),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Domestication n=514 K=2",)
dev.off()

############################ start pcadapt naive #########################

setwd("~/Documents/HG/Domestication/04_pcadapt")

# pruned SNP set with 10K as the pruning window
path_to_file <- "./genetyped_data_n_514_maf05_maxmiss07_nochr156invers.bed"
file <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./genetyped_data_n_514_maf05_maxmiss07_nochr156invers.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
par(mfrow=c(1,1))
x <- pcadapt(input = file, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("Ind514_naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
sample_file = read.delim("pop_514_sample_list.txt", header = TRUE, sep='\t')
poplist.names <- sample_file$Pop
# Computing the test statistic based on PCA
x <- pcadapt(file, K = 2)
jpeg("Ind514_naive_pca_k2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
graph2ppt(file="Ind514_naive_PC1-2.pptx", width=5, height=3)

# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
jpeg("Ind514_naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x , option = "manhattan")
dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=paste0(CHR,'_',POS), Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Ind514_naive_manhattan_by_chr.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "Domestication naive n=514 K=2",)
dev.off()

############################ start pcadapt BP #########################
setwd("~/Documents/HG/Domestication/04_pcadapt/maf005/")
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
setwd("~/Documents/HG/Domestication/04_pcadapt/maf005/")

# all
bedfile = "genetyped_data_n_514_maf05_maxmiss07_nochr156invers.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("genetyped_data_n_514_maf05_maxmiss07_nochr156invers.rds")
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
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR[outliers])
daf = data.frame(CHR=CHR, POS=POS, SNP=paste0(CHR,"_",POS), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
outlier_daf = data.frame(CHR=CHR[outliers], POS=POS[outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "Ind514_best_practice_outlier_PC1-2_10K.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_Ind514_best_practice_outlier_PC1-2_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 15),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

############################ start pcadapt BP no invers #########################
# # Full no inversions no chr56 inversions
# bedfile = "genetyped_data_n_514_maf01_maxmiss07_nochr156invers.bed"
# snp_readBed(bedfile)
# # this will create a .rds file
# obj.bigSNP <- snp_attach("genetyped_data_n_514_maf01_maxmiss07_nochr156invers.rds")
# G <- obj.bigSNP$genotypes
# SNPs <- obj.bigSNP$map$marker.ID
# CHR <- obj.bigSNP$map$chromosome
# POS <- obj.bigSNP$map$physical.pos

bedfile_noinvers = "genetyped_data_n_514_maf05_maxmiss07_nochr156invers.bed"
snp_readBed(bedfile_noinvers)
# this will create a .rds file
obj.bigSNP_noinvers <- snp_attach("genetyped_data_n_514_maf05_maxmiss07_nochr156invers.rds")
G_noinvers <- obj.bigSNP_noinvers$genotypes
SNPs_noinvers <- obj.bigSNP_noinvers$map$marker.ID
CHR_noinvers <- obj.bigSNP_noinvers$map$chromosome
POS_noinvers <- obj.bigSNP_noinvers$map$physical.pos

# obtain the "bed" snp index during pruning and manually store them into a txt file
G_noinvers <- snp_fastImputeSimple(G_noinvers, method = c("mean0"), ncores = 8)
newpc <- snp_autoSVD(G_noinvers, infos.chr = CHR_noinvers, infos.pos = POS_noinvers, thr.r2 = 0.2, size = 10) #, thr.r2 = 0.2, size = 5, min.mac = 1
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs_noinvers[which_pruned]
write.table(keep_snp_ids, file = "all_pruned_nochr156invers.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)



Gm = toMatrix(G_noinvers)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2]))
# load the na_idx
# na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
# Gm = toMatrix(G)
# Gm = Gm[,-na_idx]
# G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
# test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
print(CHR_noinvers[outliers])
daf = data.frame(CHR=CHR_noinvers, POS=POS_noinvers, SNP=paste0(CHR_noinvers,"_",POS_noinvers), Ps=-log10(result_p))
outlier_SNP = paste0(CHR_noinvers[outliers],'_',POS_noinvers[outliers])
outlier_daf = data.frame(CHR=CHR_noinvers[outliers], POS=POS_noinvers[outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
write.table(outlier_daf, file = "Ind514_best_practice_outlier_no_chr156invers_PC1-2_10K.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
jpeg("Mahattan_Ind514_best_practice_outlier_no_chr156invers_PC1-2_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 15),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

############################
###### by populaitons ######
############################
setwd("~/Documents/HG/Domestication/04_pcadapt/maf005_by_pop")
#####################
###### MEW-MES ######
#####################
toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}
pcadapt_process <- function(head){
  #head = "MEW_MES"
  bedfile_noinvers = paste0(head,".bed")
  snp_readBed(bedfile_noinvers)
  # this will create a .rds file
  obj.bigSNP_noinvers <- snp_attach(paste0(head,".rds"))
  G_noinvers <- obj.bigSNP_noinvers$genotypes
  SNPs_noinvers <- obj.bigSNP_noinvers$map$marker.ID
  CHR_noinvers <- obj.bigSNP_noinvers$map$chromosome
  POS_noinvers <- obj.bigSNP_noinvers$map$physical.pos
  
  # obtain the "bed" snp index during pruning and manually store them into a txt file
  G_noinvers <- snp_fastImputeSimple(G_noinvers, method = c("mean0"), ncores = 8)
  newpc <- snp_autoSVD(G_noinvers, infos.chr = CHR_noinvers, infos.pos = POS_noinvers, thr.r2 = 0.2, size = 10) #, thr.r2 = 0.2, size = 5, min.mac = 1
  which_pruned = attr(newpc, 'subset')
  keep_snp_ids = SNPs_noinvers[which_pruned]
  write.table(keep_snp_ids, file = paste0(head, "_all_pruned.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  Gm = toMatrix(G_noinvers)
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2]))
  # load the na_idx
  # na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
  # Gm = toMatrix(G)
  # Gm = Gm[,-na_idx]
  # G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  # test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
  result_p <- predict(test,log10 = F)
  padj <- p.adjust(result_p,method="bonferroni")
  alpha <- 0.05
  outliers <- which(padj < alpha)
  length(outliers)
  print(CHR_noinvers[outliers])
  daf = data.frame(CHR=CHR_noinvers, POS=POS_noinvers, SNP=paste0(CHR_noinvers,"_",POS_noinvers), Ps=-log10(result_p))
  outlier_SNP = paste0(CHR_noinvers[outliers],'_',POS_noinvers[outliers])
  outlier_daf = data.frame(CHR=CHR_noinvers[outliers], POS=POS_noinvers[outliers], snp_id=outlier_SNP, outlier_p = -log10(result_p)[outliers])
  write.table(outlier_daf, file = paste0(head,'_best_practice_outlier_no_chr156invers_PC1-2_10K.txt'), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  jpeg(paste0(head,"Mahattan_best_practice_outlier_no_chr156invers_PC1-2_10K.jpg"), width = 16, height = 9, units = 'in', res = 300)
  par(mfrow=c(1,1))
  manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 15),
            col=c("grey","black"),genomewideline=F, suggestiveline=F,
            ylab="-log10(p-values)", cex.lab=1.4, main = pates0(head, "PCAdapt Best Practice"))
  dev.off()
}
pcadapt_process("NCW_NCS")
