library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
library("adegenet")
library(vcfR)
library(MASS)

# see https://popgen.nescent.org/DifferentiationSNP.html for detailed example
# read from a vcf file
setwd("~/Dropbox/Mac/Documents/HG/Domestication/13_diversity_Fis_hierfstat")

# load the population information
pop_info <- read.table("pop_509_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
# load vcf file
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf --remove UMFS_2_outlier.txt --recode --recode-INFO-all --out genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe_pruned", sep=""))

#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf"
vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

# vcf_file = "genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf"
# #vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
# vcf <- read.vcfR(vcf_file, verbose = FALSE)
# Mydata2 <- vcfR2genind(vcf)
# Mydata2@pop <- pop_info$Pop_correct
# Mydata2

# Estimates allelic richness, the rarefied allelic counts, per locus and population
#Arich <- allelic.richness(Mydata1,min.n=NULL,diploid=TRUE)
#colMeans(x=Arich$Ar, na.rm = TRUE)

Arich <- allelic.richness(Mydata1,min.n=NULL,diploid=TRUE)
ind_mean <- colMeans(x=Arich$Ar, na.rm = TRUE)
# MEW1     MEW2     LIW1     LIW2     DBW1     DBW2     NCW1     NCW2     DBX1     DBX2     DBX3 
# 1.814520 1.823549 1.843482 1.844869 1.851226 1.851103 1.848713 1.821837 1.784868 1.748396 1.820937 
# UNC1     UNC2     UMFS     NEH1     NEH2     MEH2 
# 1.703802 1.689899 1.693637 1.788187 1.759862 1.646742 
wilcox.test(ind_mean[1:9],ind_mean[10:17], alternative = "greater")
mean(ind_mean[1:9])
mean(ind_mean[10:17])
# 
# Wilcoxon rank sum exact test
# 
# data:  ind_mean[1:9] and ind_mean[10:17]
# W = 69, p-value = 0.0002879
# alternative hypothesis: true location shift is greater than 0
write.table(ind_mean, file = "individual.allelic.richness.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# These statistics come from the package hierfstat Fst following Nei (1987) on genind object
basicstat <- basic.stats(Mydata1, diploid = TRUE, digits = 2) 
names(basicstat)
# $overall
# Ho   Hs   Ht  Dst  Htp Dstp  Fst Fstp  Fis Dest 
# 0.24 0.26 0.28 0.02 0.28 0.02 0.07 0.07 0.08 0.03 
# mean Ho per population
Ho <- colMeans(x=basicstat$Ho, na.rm = TRUE)
# MEW1      MEW2      LIW1      LIW2      DBW1      DBW2      NCW1      NCW2      DBX1 
# 0.2404389 0.2407640 0.2465227 0.2474728 0.2419118 0.2404444 0.2397749 0.2317948 0.2417527 
# DBX2      DBX3      UNC1      UNC2      UMFS      NEH1      NEH2      MEH2 
# 0.2373491 0.2481254 0.2457495 0.2278548 0.2306105 0.2454282 0.2454719 0.2493828 
# mean He per population
write.table(Ho, file = "pop.Ho.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

He <- colMeans(x=basicstat$Hs, na.rm = TRUE)
write.table(He, file = "pop.He.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# MEW1      MEW2      LIW1      LIW2      DBW1      DBW2      NCW1      NCW2      DBX1 
# 0.2696512 0.2706209 0.2712761 0.2721235 0.2753675 0.2760349 0.2767698 0.2673676 0.2625208 
# DBX2      DBX3      UNC1      UNC2      UMFS      NEH1      NEH2      MEH2 
# 0.2554719 0.2732786 0.2433393 0.2393306 0.2382035 0.2655503 0.2591543 0.2319194
#colMeans(x=basicstat$Fis, na.rm = TRUE)
# MEW1         MEW2         LIW1         LIW2         DBW1         DBW2         NCW1         NCW2         DBX1         DBX2         DBX3         UNC1         UNC2         UMFS 
# 0.103205353  0.107284760  0.089019027  0.085328099  0.117374068  0.122903014  0.127396284  0.125345136  0.077388928  0.072826638  0.092425782  0.007032094  0.048131693  0.047681743 
# NEH1         NEH2         MEH2 
# 0.075109389  0.054631710 -0.051158277
# Nonparametric Tests of Group Differences carried by Mann-Whitney U test
# independent 2-group Mann-Whitney U Test

wilcox.test(He[1:9],He[10:17], alternative = "greater")

# Wilcoxon rank sum exact test
# 
# data:  He[1:9] and He[10:17]
# W = 65, p-value = 0.001851
# alternative hypothesis: true location shift is greater than 0

# div <- summary(Mydata1)
# names(div)
# 
# wc(Mydata1) # Weir and Cockerham's estimate
# fst <- genet.dist(Mydata1, method = "WC84") # Pairwise Fst
# write.matrix(fst, file = "pairwise.fst.txt", ,sep = "\t")

# MEW1          MEW2          LIW1          LIW2          DBW1          DBW2          NCW1          NCW2          DBX1          DBX2          DBX3          UNC1          UNC2          UMFS          NEH1          NEH2
# MEW2 -8.338017e-04                                                                                                                                                                                                                  
# LIW1  3.397082e-02  3.041393e-02                                                                                                                                                                                                    
# LIW2  3.280786e-02  2.943828e-02  7.765745e-04                                                                                                                                                                                      
# DBW1  4.033635e-02  3.665275e-02  9.465682e-03  8.950580e-03                                                                                                                                                                        
# DBW2  4.068063e-02  3.690190e-02  9.542859e-03  9.086903e-03 -6.784757e-05                                                                                                                                                          
# NCW1  5.143801e-02  4.791863e-02  2.146300e-02  2.066784e-02  1.203849e-02  1.211993e-02                                                                                                                                            
# NCW2  6.602326e-02  6.250270e-02  3.624522e-02  3.540876e-02  2.534974e-02  2.569653e-02  3.066255e-02                                                                                                                              
# DBX1  6.666415e-02  6.262228e-02  3.695304e-02  3.640453e-02  2.756088e-02  2.741735e-02  3.892286e-02  5.286610e-02                                                                                                                
# DBX2  6.973240e-02  6.891975e-02  6.383361e-02  6.282689e-02  6.263321e-02  6.255399e-02  7.221538e-02  8.588941e-02  8.808586e-02                                                                                                  
# DBX3  4.328617e-02  4.082283e-02  4.239407e-02  4.191097e-02  4.260338e-02  4.267954e-02  4.918660e-02  6.728995e-02  6.786679e-02  6.588667e-02                                                                                    
# UNC1  1.104393e-01  1.065366e-01  8.156302e-02  8.098159e-02  7.176014e-02  7.157303e-02  6.151074e-02  8.952706e-02  9.907764e-02  1.307414e-01  1.065288e-01                                                                      
# UNC2  1.160015e-01  1.126008e-01  8.780035e-02  8.706918e-02  7.708120e-02  7.719742e-02  8.191297e-02  5.489784e-02  1.045944e-01  1.374827e-01  1.166053e-01  1.434935e-01                                                        
# UMFS  6.457170e-02  6.443535e-02  7.707049e-02  7.615761e-02  8.283863e-02  8.292444e-02  9.305922e-02  1.088684e-01  1.097428e-01  1.315982e-01  1.093249e-01  1.556493e-01  1.606666e-01                                          
# NEH1  4.698689e-02  4.353182e-02  6.444268e-02  6.365454e-02  6.924573e-02  6.911575e-02  7.574292e-02  9.293013e-02  9.418106e-02  7.904170e-02  3.664870e-02  1.335106e-01  1.426844e-01  1.297688e-01                            
# NEH2  5.104652e-02  4.721565e-02  7.072399e-02  6.976183e-02  7.573903e-02  7.581996e-02  8.433511e-02  9.970735e-02  1.010364e-01  7.295864e-02  4.674792e-02  1.426928e-01  1.498869e-01  1.342546e-01  1.380658e-02              
# MEH2  8.624712e-02  8.375016e-02  1.015456e-01  1.007302e-01  1.075827e-01  1.074138e-01  1.177968e-01  1.319653e-01  1.338656e-01  1.330166e-01  1.041087e-01  1.801793e-01  1.841926e-01  1.567893e-01  1.067608e-01  1.075763e-01

# using Kmeans and DAPC in adegenet 
# set.seed(5); dapc_a_score <- dapc(Mydata1,Mydata1$pop, n.pca = 20,n.da=10)
# temp_score <- optim.a.score(dapc_a_score)
# 
# dapc1 <- dapc(Mydata1, Mydata1$pop, n.pca = 9, n.da = 2) 
# scatter(dapc1, legend=TRUE, solid=.5) # plot of the group
# graph2ppt(file="DAPC",width=8,height=5)
# 
# percent= dapc1$eig/sum(dapc1$eig)*100
# barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))


################
#### ggplot ####
################
library(cowplot)
df <- read.delim("He_Ar_summary.csv", header = TRUE, sep=',')
df$Source <- factor(df$Source, levels=c("Wild populations", "Selected line"))
wilcox.test(df$Ne[1:8],df$Ne[9:17])

Ar <- df %>%
  ggplot(aes(Source, Ar, fill=Source)) +
  geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.1, aes(fill=Source))  +
  scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab("Allelic richness (Ar)")+
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = "none")
  #geom_point(position = "jitter", alpha = 0.7, size = 3)
He <- df %>%
  ggplot(aes(Source, He, fill=Source)) +
  geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.1, aes(fill=Source))  +
  scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab("Expected heterozygosity (He)")+
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = "none")

#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
col <- c( "#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA", 
          #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
          "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")
df$Sites = factor(df$Sites, levels=c("MEW1", "MEW2","LIW1","LIW2","DBW1","DBW2","NCW1","NCW2","DBX1","DBX2","DBX3","UNC1","UNC2","UMFS","NEH1","NEH2", "MEH2"))

Ne <- ggplot(df, aes(x=Sites, y=Ne, fill=Sites)) +
  geom_bar(stat="identity", fill=col) +
  geom_errorbar(aes(ymin=Down,ymax=Up), width=0.5)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=20),) +
  labs(x=NULL, y = "Effective population size (Ne)")  +  
  scale_fill_manual(values=col) 

N1 <- Ne + coord_cartesian(ylim = c(0, 700))
N2 <- Ne + coord_cartesian(ylim = c(10000, 50000))
up_row <- plot_grid(He, Ar)
plot_grid(
  up_row,
  N2,
  ncol = 1
)
graph2ppt(file="Diversity1",width=10,height=6)

