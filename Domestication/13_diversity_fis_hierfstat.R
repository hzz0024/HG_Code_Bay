library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
library("adegenet")
library(vcfR)

# see https://popgen.nescent.org/DifferentiationSNP.html for detailed example
# read from a vcf file
setwd("~/Dropbox/Mac/Documents/HG/Domestication/13_diversity_Fis_hierfstat")

# load the population information
pop_info <- read.table("pop_509_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "UMFS", "MEH2", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
# load vcf file
#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf"
vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_20.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf"
#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata2 <- vcfR2genind(vcf)
Mydata2@pop <- pop_info$SampleName
Mydata2

# Estimates allelic richness, the rarefied allelic counts, per locus and population
Arich <- allelic.richness(Mydata1,min.n=NULL,diploid=TRUE)
colMeans(x=Arich$Ar, na.rm = TRUE)

Arich2 <- allelic.richness(Mydata2,min.n=NULL,diploid=TRUE)
ind_mean <- colMeans(x=Arich2$Ar, na.rm = TRUE)
write.table(ind_mean, file = "individual.allelic.richness.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# MEW1     MEW2     LIW1     LIW2     DBW1     DBW2     NCW1     NCW2     UMFS     MEH2     NEH1     NEH2     DBX1     DBX2     DBX3     UNC1     UNC2 
# 1.815095 1.824171 1.843812 1.844977 1.851081 1.851103 1.848602 1.821983 1.694711 1.646872 1.787791 1.759559 1.785078 1.748407 1.820753 1.703700 1.690370 
# These statistics come from the package hierfstat Fst following Nei (1987) on genind object
basicstat <- basic.stats(Mydata1, diploid = TRUE, digits = 2) 
names(basicstat)
# $perloc
# Ho   Hs   Ht  Dst  Htp Dstp  Fst Fstp   Fis Dest
# AX.576943979 0.22 0.25 0.27 0.03 0.28 0.03 0.09 0.10  0.13 0.04
# AX.574114010 0.27 0.28 0.29 0.01 0.29 0.01 0.03 0.03  0.02 0.01
# AX.564298109 0.40 0.41 0.43 0.02 0.43 0.02 0.04 0.05  0.02 0.03
# AX.577016137 0.40 0.42 0.48 0.06 0.48 0.06 0.12 0.13  0.06 0.10
# AX.563423198 0.36 0.38 0.40 0.02 0.40 0.02 0.04 0.04  0.06 0.03
# AX.563423205 0.21 0.23 0.26 0.03 0.26 0.03 0.10 0.11  0.10 0.04
# AX.577038974 0.27 0.33 0.34 0.01 0.34 0.02 0.04 0.05  0.19 0.02
# AX.563423214 0.37 0.35 0.40 0.05 0.40 0.05 0.12 0.13 -0.05 0.08
# [ reached 'max' / getOption("max.print") -- omitted 105811 rows ]
# 
# $overall
# Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.2413 0.2617 0.2814 0.0197 0.2826 0.0210 0.0701 0.0742 0.0777 0.0284 
# mean Ho per population
colMeans(x=basicstat$Ho, na.rm = TRUE)
# MEW1      MEW2      LIW1      LIW2      DBW1      DBW2      NCW1      NCW2      UMFS      MEH2      NEH1      NEH2      DBX1      DBX2      DBX3      UNC1      UNC2 
# 0.2405095 0.2408414 0.2462831 0.2472063 0.2417180 0.2402554 0.2395589 0.2317040 0.2309214 0.2492035 0.2453556 0.2454355 0.2417120 0.2374420 0.2479245 0.2456781 0.2278786 
# mean He per population
colMeans(x=basicstat$Hs, na.rm = TRUE)
# MEW1      MEW2      LIW1      LIW2      DBW1      DBW2      NCW1      NCW2      UMFS      MEH2      NEH1      NEH2      DBX1      DBX2      DBX3      UNC1      UNC2 
# 0.2699017 0.2708632 0.2712544 0.2720264 0.2752208 0.2759058 0.2765686 0.2673267 0.2386145 0.2319808 0.2654470 0.2591021 0.2625446 0.2555047 0.2731620 0.2432748 0.2394431 
colMeans(x=basicstat$Fis, na.rm = TRUE)
# MEW1         MEW2         LIW1         LIW2         DBW1         DBW2         NCW1         NCW2         UMFS         MEH2         NEH1         NEH2         DBX1         DBX2         DBX3         UNC1         UNC2 
# 0.103423948  0.107433512  0.089160420  0.085304573  0.117272500  0.122791307  0.127165430  0.125154842  0.047818452 -0.050914240  0.075123831  0.054718835  0.077271133  0.072745055  0.092497088  0.006998047  0.048001216 



div <- summary(Mydata1)
names(div)

wc(Mydata1) # Weir and Cockerham's estimate

# $FST
# [1] 0.07107595
# 
# $FIS
# [1] 0.07917271

fst <- genet.dist(Mydata1, method = "WC84") # Pairwise Fst

# MEW1          MEW2          LIW1          LIW2          DBW1          DBW2          NCW1          NCW2          UMFS          MEH2          NEH1          NEH2          DBX1          DBX2          DBX3          UNC1
# MEW2 -8.475731e-04                                                                                                                                                                                                                  
# LIW1  3.413180e-02  3.057554e-02                                                                                                                                                                                                    
# LIW2  3.297673e-02  2.961485e-02  7.829348e-04                                                                                                                                                                                      
# DBW1  4.009651e-02  3.637702e-02  9.012688e-03  8.486938e-03                                                                                                                                                                        
# DBW2  4.048673e-02  3.667481e-02  9.106609e-03  8.641300e-03 -7.470881e-05                                                                                                                                                          
# NCW1  5.015840e-02  4.660335e-02  1.998243e-02  1.916675e-02  1.107865e-02  1.115301e-02                                                                                                                                            
# NCW2  6.354507e-02  6.006079e-02  3.361597e-02  3.280427e-02  2.340988e-02  2.374460e-02  2.832485e-02                                                                                                                              
# UMFS  6.568338e-02  6.553404e-02  7.896447e-02  7.813907e-02  8.447429e-02  8.459857e-02  9.367071e-02  1.083389e-01                                                                                                                
# MEH2  8.664105e-02  8.414885e-02  1.022093e-01  1.014512e-01  1.077700e-01  1.076438e-01  1.171004e-01  1.301621e-01  1.585445e-01                                                                                                  
# NEH1  4.724622e-02  4.379454e-02  6.501060e-02  6.425256e-02  6.947025e-02  6.938358e-02  7.513264e-02  9.112394e-02  1.315488e-01  1.070167e-01                                                                                    
# NEH2  5.136491e-02  4.754438e-02  7.137222e-02  7.044039e-02  7.601974e-02  7.614887e-02  8.370198e-02  9.791447e-02  1.360097e-01  1.078470e-01  1.378876e-02                                                                      
# DBX1  6.631960e-02  6.225187e-02  3.654619e-02  3.600623e-02  2.759592e-02  2.745506e-02  3.797121e-02  5.097833e-02  1.111488e-01  1.337566e-01  9.413359e-02  1.010194e-01                                                        
# DBX2  7.038516e-02  6.955717e-02  6.492043e-02  6.392914e-02  6.366485e-02  6.364475e-02  7.227069e-02  8.496047e-02  1.338081e-01  1.336183e-01  7.962316e-02  7.347431e-02  8.882667e-02                                          
# DBX3  4.335579e-02  4.086857e-02  4.270410e-02  4.223775e-02  4.282822e-02  4.293612e-02  4.865950e-02  6.555644e-02  1.111361e-01  1.042990e-01  3.657837e-02  4.670266e-02  6.792141e-02  6.656735e-02                            
# UNC1  1.088659e-01  1.049323e-01  8.000744e-02  7.940500e-02  7.082119e-02  7.062502e-02  6.136845e-02  8.735822e-02  1.556055e-01  1.790519e-01  1.325178e-01  1.416266e-01  9.805347e-02  1.304186e-01  1.057150e-01              
# UNC2  1.130643e-01  1.097029e-01  8.496951e-02  8.422583e-02  7.489863e-02  7.500418e-02  7.942536e-02  5.481410e-02  1.593948e-01  1.817723e-01  1.403146e-01  1.474855e-01  1.023797e-01  1.359506e-01  1.143942e-01  1.411272e-01

# using Kmeans and DAPC in adegenet 
set.seed(5); dapc_a_score <- dapc(Mydata1,Mydata1$pop, n.pca = 20,n.da=10)
temp_score <- optim.a.score(dapc_a_score)

dapc1 <- dapc(Mydata1, Mydata1$pop, n.pca = 9, n.da = 2) 
scatter(dapc1, legend=TRUE, solid=.5) # plot of the group
graph2ppt(file="DAPC",width=8,height=5)

percent= dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))

