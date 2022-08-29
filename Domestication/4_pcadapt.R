setwd("~/Dropbox/Mac/Documents/HG/Domestication/04_pcadapt")
library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(dplyr)
# to install export, one need to isntall X11(https://www.xquartz.org/), and then type devtools::install_github("tomwenseleers/export")
library(export)
library(qvalue)
library(stringr)
library(pcadapt)
library(export)
library(ggrepel)
source("manhattan.R")

#remove.packages("bigsnpr")
#install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
# load the na_idx
#na_idx <- as.integer(unlist(str_split(readLines("all_filter_idx_395.txt"), pattern = ",")))
#Gm = toMatrix(G)
#Gm = Gm[,-na_idx]
#G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
#test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))

#############################################
## Prepare dataset for PCAdapt running ######
#############################################
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

# system(paste(vcftools," --vcf DBW_NEH_n_221.recode.vcf --keep DBW_NEH_n_127.txt --recode --recode-INFO-all --out DBW_NEH_n_127", sep=""))
# system(paste(plink, " --vcf DBW_NEH_n_127.recode.vcf --allow-extra-chr --make-bed --out DBW_NEH_n_127", sep=""))
# 
# system(paste(vcftools," --vcf LIW_NEH_n_219.recode.vcf --keep LIW_NEH_n_125.txt --recode --recode-INFO-all --out LIW_NEH_n_125", sep=""))
# system(paste(plink, " --vcf LIW_NEH_n_125.recode.vcf --allow-extra-chr --make-bed --out LIW_NEH_n_125", sep=""))
# 
# system(paste(vcftools," --vcf NCW_UNC_n_104.recode.vcf --keep NCW_UNC_n_84.txt --recode --recode-INFO-all --out NCW_UNC_n_84", sep=""))
# system(paste(plink, " --vcf NCW_UNC_n_84.recode.vcf --allow-extra-chr --make-bed --out NCW_UNC_n_84", sep=""))
# 
# system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep DBW_LIW_NEH_NCW_UNC_n_300.txt --recode --recode-INFO-all --out DBW_LIW_NEH_NCW_UNC_n_300", sep=""))
# system(paste(plink, " --vcf DBW_LIW_NEH_NCW_UNC_n_300.recode.vcf --allow-extra-chr --make-bed --out DBW_LIW_NEH_NCW_UNC_n_300", sep=""))
# 
# system(paste(vcftools," --vcf NCW_UNC_n_104.recode.vcf --keep NCW_UNC_n_84.txt --recode --recode-INFO-all --out NCW_UNC_n_84", sep=""))
# system(paste(plink, " --vcf NCW_UNC_n_84.recode.vcf --allow-extra-chr --make-bed --out NCW_UNC_n_84", sep=""))
# 
# system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep DBW_LIW_NEH_NCW_UNC_n_300.txt --recode --recode-INFO-all --out DBW_LIW_NEH_NCW_UNC_n_300", sep=""))
# system(paste(plink, " --vcf DBW_LIW_NEH_NCW_UNC_n_300.recode.vcf --allow-extra-chr --make-bed --out DBW_LIW_NEH_NCW_UNC_n_300", sep=""))

#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep no_MEH_MEW_n_416.txt --maf 0.05 --recode --recode-INFO-all --out no_MEH_MEW_n_416", sep=""))
#system(paste(plink, " --vcf no_MEH_MEW_n_416.recode.vcf --allow-extra-chr --make-bed --out no_MEH_MEW_n_416", sep=""))

#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep no_DBX1_UNC_MEH_MEW_MEH_n_342.txt --maf 0.05 --recode --recode-INFO-all --out no_DBX1_UNC_MEH_MEW_MEH_n_342", sep=""))
#system(paste(plink, " --vcf no_DBX1_UNC_MEH_MEW_MEH_n_342.recode.vcf --allow-extra-chr --make-bed --out no_DBX1_UNC_MEH_MEW_MEH_n_342", sep=""))

system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_477.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_477", sep=""))
system(paste(plink, " --vcf pop_n_477.recode.vcf --allow-extra-chr --make-bed --out pop_n_477", sep=""))

system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_282.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_282", sep=""))
system(paste(plink, " --vcf pop_n_282.recode.vcf --allow-extra-chr --make-bed --out pop_n_282", sep=""))

system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_123.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_123", sep=""))
system(paste(plink, " --vcf pop_n_123.recode.vcf --allow-extra-chr --make-bed --out pop_n_123", sep=""))


system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_104.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_104", sep=""))
system(paste(plink, " --vcf pop_n_104.recode.vcf --allow-extra-chr --make-bed --out pop_n_104", sep=""))


#############################################
## Function for PCadapt best practise (BP) ##
#############################################

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}
# colours   DBW1      DBW2         DBX1      DBX2         DBX3      LIW1        LIW2        MEH2       MEW1      MEW2         NCW1      NCW2         NEH1      NEH2       UMFS       UNC1       UNC2
col <- c("#33A65C", "#57A337",  "#F64971", "#FC719E", "#EB73B3",  "#30BCAD", "#21B087",  "#F89217", "#1BA3C6", "#2CB5C0",  "#A2B627", "#F8B620", "#F06719", "#E03426",  "#F8B620", "#CE69BE", "#A26DC2")

sub_name = "no_MEH_MEW_n_416"
# colours   DBW1      DBW2         DBX1      DBX2         DBX3      LIW1        LIW2       NCW1      NCW2         NEH1      NEH2       UMFS       UNC1       UNC2
col <- c("#33A65C", "#57A337",  "#F64971", "#FC719E", "#EB73B3",  "#30BCAD", "#21B087",  "#A2B627", "#F8B620", "#F06719", "#E03426",  "#F8B620", "#CE69BE", "#A26DC2")
alpha = 0.05

sub_name = "no_DBX1_UNC_MEH_MEW_MEH_n_342"
# colours   DBW1      DBW2         DBX2       DBX3       LIW1        LIW2        NCW1       NCW2         NEH1      NEH2       UMFS       
col <- c("#ec9c9d", "#8ea4bf")
#col <- c("#87CEEB", "#1E90FF",  "#FC719E", "#EB73B3",  "#778899", "#C0C0C0",  "#98FB98", "#2E8B57", "#F06719", "#E03426",  "#F8B620")
alpha = 0.05

sub_name = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

sub_name = "pop_n_282"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

sub_name = "pop_n_104"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

sub_name = "pop_n_123"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

sub_name = "pop_n_477"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

pcadapt_BP <- function(sub_name, alpha){
  
  f_bk = paste0(sub_name, ".bk")
  if (file.exists(f_bk)) {
    #Delete file if it exists
    file.remove(f_bk)
  }
  
  # part 1 SNP clumpping and data preparation
  snp_readBed(paste0(sub_name, ".bed"))
  # this will create a .rds file
  obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
  G <- obj.bigSNP$genotypes
  SNPs <- obj.bigSNP$map$marker.ID
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  # check if there is any missing values as NA
  #big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
  # genotype imputation
  G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
  #big_counts(G, ind.col = 1:dim(G)[1]) # check if NAs are 0
  # LD clumping using r2 = 0.2
  newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) # size is the window size of 10K
  # extract SNPs after clumpping
  which_pruned = attr(newpc, 'subset')
  keep_snp_ids = SNPs[which_pruned]
  write.table(keep_snp_ids, file = paste0(sub_name, "_pruned_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  print(paste0("SNPs after clumpping is: ", length(keep_snp_ids), " out of ", dim(obj.bigSNP$map)[1]))
  # generate thinned vcf file
  system(paste(vcftools," --vcf ",sub_name,".recode.vcf", " --snps ", sub_name, "_pruned_SNP.txt", " --remove remove.txt --recode --recode-INFO-all --out ", sub_name, "_pruned", sep=""))
  system(paste(plink, " --vcf ",sub_name,"_pruned.recode.vcf", " --allow-extra-chr --make-bed --out ", sub_name,"_pruned", sep=""))
  
  # part 2 Population structure patterns from thinned data
  file <- read.pcadapt(paste0(sub_name,"_pruned.bed"), type = "bed")
  bim_file = read.delim(paste0(sub_name,"_pruned.bim"), header = FALSE, sep='\t')
  # Scree plot
  x <- pcadapt(input = file, K = 10)
  # The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
  scree_plot <- plot(x, option = "screeplot")
  # Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
  sample_file = read.delim(paste0(sub_name,"_list.txt"), header = TRUE, sep='\t')
  poplist.names <- sample_file$Pop
  # PCA plot
  PCA_plot <- plot(x, option = "scores", i=1, j=2, pop=poplist.names, col=col)
  
  # part 3 Best practise for outlier detection
  Gm = toMatrix(G)
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  Gm = Gm[,snp_MAF(G_coded) != 0] # eliminate monomorphic SNPs
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  #res_pcadapt <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2])) # use PC1 and PC2 for outlier identification
  res_pcadapt <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2])) # use defined PCs for outlier identification
  res_p <- predict(res_pcadapt,log10 = F)
  hist(res_p, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
  # Choosing a cutoff for outlier detection
  qval <- qvalue(res_p)$qvalues
  outliers <- which(qval < alpha)
  print(paste0("number of outliers from best practise in ", sub_name, " is: ", length(outliers)))
  # start to make Mahattan plot
  CHR = CHR[which(snp_MAF(G_coded) != 0)]
  POS = POS[which(snp_MAF(G_coded) != 0)]
  daf = data.frame(CHR=CHR, POS=POS, SNP=paste0(CHR,"_",POS), Ps=res_p, qvalue=qval)
  outlier_SNP <- paste0(CHR[outliers],'_',POS[outliers])
  outlier_daf <- data.frame(CHR=CHR[outliers], POS=POS[outliers], POS=POS[outliers], snp_id=outlier_SNP)
  write.table(outlier_daf, file = paste0(sub_name, "_PCAdapt_BP_q05_n_", length(outliers), ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  #jpeg("Mahattan_Ind509_best_practice_outlier_PC1-2_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
  # customize mahattan plot see https://www.r-graph-gallery.com/101_Manhattan_plot.html
  don <- daf %>% 
    dplyr::group_by(CHR) %>% 
    dplyr::summarise(chr_len=max(POS)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>% 
    left_join(daf, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, POS) %>%
    mutate(BPcum=POS+tot) %>%
    mutate(is_highlight=ifelse(SNP %in% outlier_SNP, "yes", "no"))
  axisdf = don %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  write.table(don, file = paste0(sub_name, "_qs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  outlier_plot <- ggplot(don, aes(x=BPcum, y=-log10(qvalue))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
    scale_color_manual(values = rep(c("grey", "#4F4F4F"), 22 )) +
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0.2)) + # remove space between plot area and x axis
    geom_point(data=subset(don, is_highlight=="yes"), color="red", size=0.8) +
    xlab("Chromosome") + 
    ylab("-log10(qvalue)")+
    theme_bw() +
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    ylim(0, 9) +  # for better visualization, need to remove this y scale limit for complete plot
    labs(caption=paste0(dim(don[-log10(don$qvalue) > 15,])[1], " out of ", length(outliers), " outliers not shown due to Y-axis setting < 15"))
  
  bottom_row <- plot_grid(scree_plot, PCA_plot, labels = c('A', 'B'), align = 'h', rel_widths = c(1, 1.3))
  final_plot <- plot_grid(bottom_row, outlier_plot, labels = c('', 'C'), ncol = 1, rel_heights = c(1, 1.2))
  ggsave(paste0("figs/", sub_name, ".jpg"), final_plot, width = 9, height = 7)
}

pcadapt_BP("genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe", 0.05)
#pcadapt_BP("no_DBX1_UNC_MEH_MEW_MEH_n_342", 0.05)
pcadapt_BP("pop_n_282", 0.05)
pcadapt_BP("pop_n_104", 0.05)
pcadapt_BP("pop_n_123", 0.05)

pcadapt_BP("pop_n_477", 0.05)
#pcadapt_BP("no_MEH_MEW_n_416", 0.05)

######## test ########

MY2 <- pcadapt::read.pcadapt("DBW_NEH_n_127.recode.vcf", type = "vcf")
sample_file = read.delim("DBW_NEH_n_127_list.txt", header = TRUE, sep='\t')
poplist.names <- sample_file$Pop
x <- pcadapt(MY2,K=20)
x2<-pcadapt(MY2,K=1)
plot(x,option="screeplot")
plot(x, option = "scores", pop = poplist.names)
plot(x, option="manhattan")
plot(x, option = "qqplot", threshold = 0.1)
hist(x2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x2, option = "stat.distribution")
qval <- qvalue(x2$pvalues)$qvalues
outliers2 <- which(qval<alpha)
length(outliers2)
mat <- mat <- mat <- pcadapt::bed2matrix(MY2)
hist(colMeans(mat, na.rm = TRUE))
hist(colMeans(is.na(mat)))
hist(rowMeans(is.na(mat)))
N <- nrow(mat)
M <- 50e3
mat_null <- sapply(runif(M, min = 0.05, max = 0.5), function(af) {
  rbinom(N, size = 2, prob = af)
})
mat2 <- cbind(mat, mat_null)
MY2 <- read.pcadapt(mat2, type = "lfmm")


system("bedtools intersect -a test1.bed -b test2.bed > test.intersect", intern = TRUE)


