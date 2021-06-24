# I did not explore other functions of DAPC (e.g. check allele changes, group assignment)

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
# load the input (genepop format with the extension of .gen)
# set your working directory
setwd("load your working directory")
# install the package "adegenet"
install.packages("adegenet")
# load adengnet
library(adegenet)
# this will open DAPC server and allow you to initially visualize your data structure
# play with this server to find the approximate cluster numbers in your data
adegenetServer(what = "DAPC")
# reading the genepop file as input
x <- read.genepop("DP3g98maf0105cleanHW.SNP2_10dist.chr6.2pop.gen",quiet = TRUE)
x
# print and check samples and poopulations
x@pop

###############################################
# generate initial DAPC data
dapc1 <- dapc(x, x@pop)
# when asked "Choose the number PCs to retain (>=1):"
# I would choose first few PCs that explain the majortiy of cumulated variance, here is 10
# when asked "Choose the number discriminant functions to retain (>=1):"
# I would choose first few DAs that have high F-statistic values, here is 2

###############################################
# identifying genetic clusters using find.clusters
grp <- find.clusters(x, max.n.clust=25)
# when asked "Choose the number PCs to retain (>=1):"
# I would choose first few PCs that explain the majortiy of cumulated variance, here is 10
# when asked "Choose the number of clusters (>=2:"
# here the function displays a graph of BIC values for increasing values of k
# This graph shows a significant decrease of BIC until k = 4 clusters.
# In this case, the elbow in the curve also matches the smallest BIC, clearly indicates 4 clusters should be retained.

###############################################
# Alpha score test for optimal number of PCs
# option test for alpha score with 20 simulations, this function allow you examine what is the optiomal number of PCs to retain
temp <- optim.a.score(dapc1, n.sim = 20)

# load a package called "export" for powerpoint file export (an easy way to edit or combine your plots)
library(export)
graph2ppt(file="WAL605_SNP68_alpha",width=4,height=6)

# generate DAPC data again with the best optimal number of PCs from alpha-score evaluation
# or you can just use whatever PCs you want by specify the n.pca = 10 (i.e. using first 10 PCs)
dapc1 <- dapc(x, x@pop, n.pca = temp$best)
# when asked "Choose the number discriminant functions to retain (>=1):"
# I would choose first few DAs that have high F-statistic values, here is 2

###############################################
# DAPC plot

# define color class, here is rainbow
col1 <- rainbow(length(levels(pop(x))))

# within the scatter function, you can
# edit cstar=0 will remove the connection lines
# edit cell=0  will remove the ellipses
# edit cex=1.5 to adjust the size of sample dots
# add label.inds = list(air = 2, pch = NA) to label the dots
scatter(dapc1,  xlab="Discriminatn Function 1", ylab="Discriminatn Function", cex = 1.5, pch=15:18, lwd=2, lty=2, solid = .8, cstar = 1, cell = 1, scree.da=FALSE,  col = col1, clab = 0 )
# plot the legend
legend("topleft", legend = levels(pop(x)),cex = 0.5,col = col1, bty = 'n', pch=15:18, inset = c(0,0), xpd = TRUE)
# plot the PCA variance 
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=15, type="h", lwd=3)
}
add.scatter(myInset(), posi="bottomleft",
            inset=c(0.1,0), ratio=.14,
            bg=transp("white"))

# again load a package called "export" for powerpoint file export (an easy way to edit or combine your plots)
library(export)
graph2ppt(file="WAL605_SNP68_DAPC",width=8,height=6)


