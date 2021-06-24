# import genepop to adegenet and make basic PCA plot
# need code to calculate AR/AR, AR/ST, ST/ST inversion genotype frequencies
# do Hardy-Weinberg test on genotype frequencies
install.packages("wordcloud")
###############################################
# load the input (genepop format with the extension of .gen)
# set your working directory
setwd("load your working directory")
# install the package "adegenet"
#install.packages("adegenet")
# load adengnet
library(adegenet)
# this will open DAPC server and allow you to initially visualize your data structure
# play with this server to find the approximate cluster numbers in your data
#adegenetServer(what = "DAPC")
# reading the genepop file as input
obj1 <- read.genepop("DP3g98maf0105cleanHW.SNP2_10dist.chr6.gen",quiet = TRUE, ncode=2) #ncode is number of digits in allele designation
obj1
# print and check samples and poopulations
obj1@pop

###############################################

#Allele presence absence data are extracted and NAs replaced using tab:
X <- tab(obj1, NA.method="mean")
X

## make PCA
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(obj1))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
## basic plot
plot(pca1$li, col=myCol, cex=1, pch=myPch, main = "DP3g98maf0105cleanHW.SNP2_10dist.chr6")

## write file of per-individual loadings
write.csv(pca1$li, file = "DP3g98maf0105cleanHW.SNP2_10dist.chr6.loadings.csv")
## use wordcloud for non-overlapping labels
library(wordcloud)
## Loading required package: RColorBrewer
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=0.5, new=TRUE)
## legend the axes by adding loadings
abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1*.5, add.plot=TRUE)
legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),
       leg=c("Hastings","Soundview"), pt.cex=2)
###########################################################################
## chi-square test, first read the PC-axis table 
dt <- read.csv(file = "DP3g98maf0105cleanHW.SNP2_10dist.chr6.loadings.csv")
## extract the positions of first and last HH samples, change the "^H" to "^S" for SV population
first <- grep("^H", dt[,1])[1]
last <- tail(grep("^H", dt[,1]), n=1)
## count the observed ST, HE(AR x ST), and AR, here -5 and 1.2 are mannually assigned for three cluster distinction
ST <- sum((dt[first:last,2]) < -5) # change the -5 here for different PCA plot
HE <- sum((dt[first:last,2]) < 1.2 & (dt[first:last,2]) > -5) # change the 1.2 (x-axis position of HH912s_121) and -5 here for different PCA plot
AR <- sum((dt[first:last,2]) > 1.2) # change the 1.2 here for different PCA plot
## total observation counts
total <- sum(ST,HE,AR)
## ST and AR frequency calculation
ST_freq <- ((ST+(0.5*HE))/total)
AR_freq <- ((AR+(0.5*HE))/total)
## expected ST, HE(AR x ST), and AR counts under HWE
exp_ST <- ST_freq^2*total
exp_AR <- AR_freq^2*total
exp_HE <- 2*ST_freq*AR_freq*total
## chi-square test, the result is p-value
Chi_Sq = ((ST-exp_ST)^2 / exp_ST) + ((HE-exp_HE)^2 / exp_HE) + ((AR-exp_AR)^2 / exp_AR)
pchisq(Chi_Sq, df=1, lower.tail=FALSE)

