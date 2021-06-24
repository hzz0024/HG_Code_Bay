#################################################
# Script for computing the neutral distribution #
# of "delta p" values under neutrality and      #
# finite sampling with missing genotypes        #
#################################################


## for N observed diploid genotypes, X is the number of counts ##
## of the reference allele in [1,2N-1].                        ##
## For a given allele frequency p in the gene pool, the proba  ##
## of observing X alleles in N diploid individuals is given by ##  
## the density of the binomial distribution:                   ##
## P(X|p) = dbinom(X,2N,p)                                     ##
## Now we want the probability density of p in the gene pool   ##
## given X observed alleles in N diploid individuals. This is  ##
## obtained using Bayes theorem which gives:                   ##
## P(p|X) = (P(X|p)*P(p))/P(X)                                 ##
## in which P(p) is obtained using equation 50 in Tajima 1989  ##

#######################
# ESTIMATION OF THETA #
#######################

## We use a function that estimates theta from our data using  ##
## the mean number of SNPs (s) per polymorphic RAD locus after ##
## quality filtering.                                          ##

### PARAMETERS ###

# -l the length of RAD tags in bp
# -n the total number of individuals used to calculate -s
# -s the mean number of SNPs per RAD

Theta <- function(s,l,n)
  {
  # from equation (5) in Tajima 1989 we have M = S/ai, with ai defined as
  ai = NULL
  for (i in 1:((2*n)-1))
    {ai = c(ai,1/i)}
  ai = sum(ai)
  # So using the number of SNPs per pb S = (s/l) we get theta as
  M = (s/l)/ai
  }

theta <- Theta(2.185,86,251)

## In our case this gives theta = 0.00374 which is close to Pi ##

######################
# AFS AT EQUILIBRIUM #
######################

## Now we use equation 50 in Tajima 1989 to get the Allele     ##
## Frequency Spectrum (Pp) over [1,2N-1], which sum is scaled  ##
## to one to get the probability density of allele counts      ##

### PARAMETERS ###

# -theta the value of theta estimated using the function Theta()
# -N the number of individuals that were sampled 

AFS <- function(theta,N)
  {
  Pp = NULL
  for (j in 1:(2*N-1))
    {Pp = c(Pp,theta*((1/j)+(1/(2*N-j))))}
  Pp = Pp*(1/sum(Pp))
  }


#######################################################
# PROBABILITY DENSITY OF THE UNKNOWN ALLELE FREQUENCY #
#######################################################

## So now we come back to P(p|X) = (P(X|p)*P(p))/P(X)  ##

### PARAMETERS ###

# -theta the value of theta estimated using the function Theta()
# -N the number of individuals that were sampled 
# -X the count of the reference allele in the sample of size N

Dfreq <- function(X,N,theta)
  {
  # p is the vector that stores all the frequencies that the
  # allele can take in a sample of size N
  p=(1:(2*N-1))/(2*N)
  # we get P(X|p) over all possible values of p within the pool
  PXp = dbinom(X,2*N,p)
  # we get the probability density of allele counts P(p):
  Pp = AFS(theta,N)
  # and P(X) over [1,2N-1] is simply P(X) = 1/(2N-1)
  PX = 1/(2*N-1)
  # So for a given value of X observed in a sample of size N
  # we get P(p|X) the probability density of p in the gene pool
  PpX = NULL
  for(k in 1:(2*N-1))
    {PpX=c(PpX,(PXp[k]*Pp[k])/PX)}
  # and we scale it to one
  PpX = PpX*(1/sum(PpX))
  }

plot(Dfreq(25,50,theta))  


################################
# PROBABILITY DENSITY OF DELTA #
################################

## Now we randomly draw allele frequency values from ##
## the probability density of the unknown allele     ##
## frequency in the gene pool to make samples of     ##
## size K using binomial sampling                    ##

### PARAMETERS ###

# -theta the value of theta estimated using the function Theta()
# -N the number of individuals that were sampled at T0 before selection
# -X the count of the reference allele in the sample of size N
# -K the number of individuals that were sampled at T1 after selection
# -nreps the number of random sampling events to draw the neutral distribution of delta

Ddelta <- function(X, N, K, theta, nreps)
  {
  poolfreq = colSums((1:(2*N-1))*rmultinom(nreps,1:(2*N-1),prob=matrix(Dfreq(X,N,theta),1)))/(2*N)
  samplefreq = rbinom(nreps,2*K,poolfreq)/(2*K)
  delta = (X/(2*N)) - samplefreq
  }


##########################
# DELTA P LEVENE'S MODEL #
##########################

## We calculate the value of delta P after selection   ##
## in one habitat under Levene's moldel, as a function ##
## of initial allele frequency before selection and    ##
## the strength of selection. The fitness of AA is 1+s ##
## the fitness of Aa is 1 and the fitness of aa is 1-s ##
## Therefore, allele A determines genotypes' fitness   ##
## additively, and selection has an antagonistic and   ##
## symmetrical effect on the fitness of homozygotes    ##


### PARAMETERS ###

# (-P internal parameter) the vector of allele frequencies in the habitat (i.e. in the larval pool) before selection
# -s the selection coefficient s

Levene <- function(s)
  {
  P <- seq(0,1,0.01)
  #DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s)
  DpLev <- s*P*(1-P)/(1 + 2*s*P - s) #
  }


library(ggplot2)

p <- rep(seq(0,1,0.01),101)
s <- sort(rep(seq(0,0.5,0.005),101))
delta <- unlist(lapply(seq(0,0.5,0.005),Levene))
data <- data.frame(p,s,delta)
d <- ggplot(data, aes(p, s, z = delta))
d + stat_summary2d() + scale_fill_gradient2(name = "Delta p", high = "red", mid = "grey",low = "blue",midpoint = 0.05,guide = "colourbar")


####### change by HG
#p <- rep(seq(0,1,0.01),101)
p <- rep(seq(0,0.99,0.01),100)
s <- sort(rep(seq(0,0.99,0.01),100))
#s <- sort(rep(seq(0,0.495,0.005),100)) 
#delta <- unlist(lapply(seq(0,0.495,0.005),Levene))
delta <- unlist(lapply(seq(0,0.99,0.01),Levene))
data <- data.frame(p,s,delta)
write.table(data, "./data.csv", sep=",", quote=F, row.names=T, col.names=NA)
d <- ggplot(data, aes(p, s, z = delta))
d + stat_summary2d() + scale_fill_gradient2(name = "Delta p", high = "red", mid = "grey",low = "black",midpoint = 0.1,guide = "colourbar")



##############################
# POWER OF THE BAYESIAN TEST #
##############################

## We use the initial and post-selection frequencies ##
## calculated with the Levene function and stored in ##
## data. First, the initial allele frequency is used ##
## to build the pre-selection sample of size N. Then ##
## a random perturbation is applied to the allele    ##
## frequencies post selection to model the effect of ##
## drift in a finite population of size 10000. Then  ##
## these new post selection frequencies are used to  ##
## build a post-selection sample of size K.          ##

### PARAMETERS ###

# -N the number of individuals that were sampled at T0 before selection
# -K the number of individuals that were sampled at T1 after selection
# -data the dataframe build with the Levene function


Power <- function(N, K, data)
  {
  d <- data
  #d$power <- NA
  d$powerfish <- NA
  for (r in (1:nrow(data)))
    {
    IAC <- rbinom(100,(2*N),data[r,1])
    IAF <- IAC/(2*N)
    SAC <- rbinom(100,(2*K),(rbinom(1,20000,(data[r,1]+data[r,3]))/20000))
    SAF <- SAC/(2*K)
    DAF <- IAF-SAF
    #PV <- NULL
    PVfish <- NULL
    for (t in (1:100))
      {
      if (DAF[t]==0)
        {
        #PV <- c(PV, 1)
        PVfish <- c(PVfish, 1)
        }
      if (DAF[t]>0)
        {
        #PV <- c(PV, length(which(Ddelta(IAC[t],N,K,theta,10000)>=DAF[t]))/10000)
        PVfish <- c(PVfish, fisher.test(matrix(c(IAC[t],(2*N)-IAC[t],SAC[t],(2*K)-SAC[t]),nrow=2))$p)
        }
      if (DAF[t]<0)
        {
        #PV <- c(PV, length(which(Ddelta(IAC[t],N,K,theta,10000)<=DAF[t]))/10000)
        PVfish <- c(PVfish, fisher.test(matrix(c(IAC[t],(2*N)-IAC[t],SAC[t],(2*K)-SAC[t]),nrow=2))$p)
        }
      }
    #d$power[r] <- length(which(PV<0.05))
    d$powerfish[r] <- length(which(PVfish<0.05))
    print(r)
    }
  d
  }

POW <- Power(40,100,data)

POW <- Power(50,50,data)
write.table(POW, "./POW50_orig_fish_s0_1.csv", sep=",", quote=T, row.names=T, col.names=NA)

# Summary plot (N0=40,N1=100) #
D <- ggplot(POW, aes(p, s, z = power))
D <- D + stat_summary2d() + theme_bw() + labs(x="",y="")
D + scale_fill_gradient2(name="",limits=c(0,80),breaks=c(0,20,40,60,80),high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")

E <- ggplot(POW, aes(p, s, z = powerfish))
E <- E + stat_summary2d() + theme_bw() + labs(x="",y="")
E + scale_fill_gradient2(name="",limits=c(0,80),breaks=c(0,20,40,60,80),high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")

G <- ggplot(POW, aes(p, s, z = (power-powerfish)))
G <- G + stat_summary2d() + theme_bw() + labs(x="",y="")
G +  scale_fill_gradient2(name="",limits=c(0,20),breaks=c(0,5,10,15,20),high="red",mid="grey",low="blue",midpoint=10,guide="colourbar")

###################################### formal plot ###################################### 
install.packages("ggplot2")
library("ggplot2")
library("gridExtra")
require(grid)

file2 = "POW50_orig_fish_s0_1.csv"
POW50 <- read.delim(file2, header = TRUE, sep=',')
max(POW50$delta)
#Levene delta p
d <- ggplot(data, aes(p, s, z = delta))
d <- d + stat_summary2d() + theme_bw() + labs(x="p",y="s")
deltap <- d + stat_summary2d() + scale_fill_gradient2(name = expression(Delta~italic(p)), limits=c(0,0.5), high = "red", mid = "grey",low = "#B6E5D8",midpoint = 0.25,guide = "colourbar") +
          theme(text = element_text(size=20))
deltap
# data delta p
D2 <- ggplot(POW50, aes(p, s, z=delta))
D2 <- D2 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_0 <- D2+ scale_fill_gradient2(name=expression(Delta~italic(p)),limits=c(0,0.45),breaks=c(0,0.1,0.2,0.3, 0.4),high="red",mid="grey",low="#B6E5D8",midpoint=0.2,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_0

D2 <- ggplot(POW50, aes(p, s, z = powerfish))
D2 <- D2 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_1 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_1

D2 <- ggplot(POW50, aes(s , delta , z = powerfish))
D2 <- D2 + stat_summary2d() + theme_classic() + labs(x="Selection coefficient",y=expression(Delta~italic(p)))
plot50_2 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_2

D2 <- ggplot(POW50, aes(p , delta , z = powerfish))
D2 <- D2 + stat_summary2d() + theme_classic() + labs(x=expression(italic(p)*0),y=expression(Delta~italic(p)))
plot50_3 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_3

grid.arrange(plot50_0, plot50_1, plot50_2, plot50_3, nrow = 2)
grid_arrange_shared_legend(plot50_2, plot50_3, position = "right")

# function for merging legend
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }














POW100 <- Power(100,100,data)

# Summary plot (N0=100,N1=100) #
D <- ggplot(POW100, aes(p, s, z = power))
D <- D + stat_summary2d() + theme_bw() + labs(x="",y="")
D + scale_fill_gradient2(name="",limits=c(0,100),breaks=c(0,25,50,75,100),high="red",mid="grey",low="blue",midpoint=50,guide="colourbar")

E <- ggplot(POW100, aes(p, s, z = powerfish))
E <- E + stat_summary2d() + theme_bw() + labs(x="",y="")
E + scale_fill_gradient2(name="",limits=c(0,100),breaks=c(0,25,50,75,100),high="red",mid="grey",low="blue",midpoint=50,guide="colourbar")

G <- ggplot(POW100, aes(p, s, z = (power-powerfish)))
G <- G + stat_summary2d() + theme_bw() + labs(x="",y="")
G +  scale_fill_gradient2(name="",limits=c(0,20),breaks=c(0,5,10,15,20),high="red",mid="grey",low="blue",midpoint=10,guide="colourbar")


# detailed plots #
D <- ggplot(POW, aes(p, s, z = power))
D + stat_summary2d(bins=101) + theme_bw() + scale_fill_gradient2(name="Power",high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")
E <- ggplot(POW, aes(p, s, z = powerfish))
E + stat_summary2d(bins=101) + theme_bw() + scale_fill_gradient2(name="Power Fisher",high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")
G <- ggplot(POW, aes(p, s, z = (power-powerfish)))
G + stat_summary2d(bins=101) + theme_bw() + scale_fill_gradient2(name="Delta power",high="red",mid="grey",low="blue",midpoint=10,guide="colourbar")


# setting the working directory
setwd ("C:/Users/Pierre-Alexandre/Documents/Recherche/Daurades/Delta_freq")
# load SNP dataset
SNP <- read.table("batch_2_couv_hwe_ok_snp.txt",header=TRUE)

RES <- NULL
for (u in 2:ncol(SNP))
  {
  AAL <- length(which(SNP[208:251,u]==0))
  ABL <- length(which(SNP[208:251,u]==1))
  BBL <- length(which(SNP[208:251,u]==2))
  AFL <- (2*AAL+ABL)/(2*(AAL+ABL+BBL))
  AAM <- length(which(SNP[1:105,u]==0))
  ABM <- length(which(SNP[1:105,u]==1))
  BBM <- length(which(SNP[1:105,u]==2))
  AFM <- (2*AAM+ABM)/(2*(AAM+ABM+BBM))
  AAT <- length(which(SNP[106:207,u]==0))
  ABT <- length(which(SNP[106:207,u]==1))
  BBT <- length(which(SNP[106:207,u]==2))
  AFT <- (2*AAT+ABT)/(2*(AAT+ABT+BBT))
  # Test delta Larves/Mauguio
  if (AFL-AFM==0)
  {PVM <- 1}
  if (AFL-AFM>0)
  {PVM <- length(which(Ddelta((2*AAL+ABL),(AAL+ABL+BBL),(AAM+ABM+BBM),theta,10000)>=(AFL-AFM)))/10000}
  if (AFL-AFM<0)
  {PVM <- length(which(Ddelta((2*AAL+ABL),(AAL+ABL+BBL),(AAM+ABM+BBM),theta,10000)<=(AFL-AFM)))/10000}
  # Test delta Larves/Thau
  if (AFL-AFT==0)
  {PVT <- 1}
  if (AFL-AFT>0)
  {PVT <- length(which(Ddelta((2*AAL+ABL),(AAL+ABL+BBL),(AAT+ABT+BBT),theta,10000)>=(AFL-AFT)))/10000}
  if (AFL-AFT<0)
  {PVT <- length(which(Ddelta((2*AAL+ABL),(AAL+ABL+BBL),(AAT+ABT+BBT),theta,10000)<=(AFL-AFT)))/10000}
  # Test delta Mauguio/Thau
  if (AFM-AFT==0)
  {PVJ <- 1}
  if (AFM-AFT>0)
  {PVJ <- length(which(Ddelta((2*AAM+ABM),(AAM+ABM+BBM),(AAT+ABT+BBT),theta,10000)>=(AFM-AFT)))/10000}
  if (AFM-AFT<0)
  {PVJ <- length(which(Ddelta((2*AAM+ABM),(AAM+ABM+BBM),(AAT+ABT+BBT),theta,10000)<=(AFM-AFT)))/10000}
  # AFL 90% CI
  AFL90 <- c(max((1:(2*(AAL+ABL+BBL)-1))[cumsum(Dfreq((2*AAL+ABL),(AAL+ABL+BBL),theta))<0.05])/(2*(AAL+ABL+BBL)),min((1:(2*(AAL+ABL+BBL)-1))[cumsum(Dfreq((2*AAL+ABL),(AAL+ABL+BBL),theta))>0.95])/(2*(AAL+ABL+BBL)))
  if (AFL==1) {AFL90[2]=1}
  # Store stats
  RES <- rbind(RES, c((AAL+ABL+BBL),(AAM+ABM+BBM),(AAT+ABT+BBT),AFL,AFL90,AFM,AFT,(AFL-AFM),(AFL-AFT),(AFM-AFT),PVM,PVT,PVJ))
  print(u)
  }

write.table(RES,"SGS_test_out.txt",sep="\t")


OUT <- read.table("Outliers_0.001_SGS.txt",header=TRUE)
plot(OUT$p_Larves,col="red",pch="-",cex=1.5,ylim=c(0.4,1),ylab=expression(italic("p")),xlab="Loci",xlim=c(0,70))
for (l in 1:nrow(OUT))
  {
  segments(l,OUT$pL_CI_low[l],l,OUT$pL_CI_high[l],col="grey",lwd=1.5)
  points(OUT$p_Larves,col="red",pch="-",cex=1.5)
  if (OUT$Delta_p_M[l]<=0)
    {
    if (OUT$PV_M[l]>0.05)
      {points(l,OUT$p_Mauguio[l],col="green",pch=2,cex=1.2)}
    if (OUT$PV_M[l]<=0.05)
      {points(l,OUT$p_Mauguio[l],col="green",pch=24,bg="green",cex=1.1)}
    }
  if (OUT$Delta_p_T[l]>=0)
    {
    if (OUT$PV_T[l]>0.05)
      {points(l,OUT$p_Thau[l],col="blue",pch=6,cex=1.2)}
    if (OUT$PV_T[l]<=0.05)
      {points(l,OUT$p_Thau[l],col="blue",pch=25,bg="blue",cex=1.1)}
    } 
  if (OUT$Delta_p_M[l]>0)
    {
    if (OUT$PV_M[l]>0.05)
      {points(l,OUT$p_Mauguio[l],col="green",pch=6,cex=1.2)}
    if (OUT$PV_M[l]<=0.05)
      {points(l,OUT$p_Mauguio[l],col="green",pch=25,bg="green",cex=1.1)}
    }
  if (OUT$Delta_p_T[l]<0)
    {
    if (OUT$PV_T[l]>0.05)
      {points(l,OUT$p_Thau[l],col="blue",pch=2,cex=1.2)}
    if (OUT$PV_T[l]<=0.05)
      {points(l,OUT$p_Thau[l],col="blue",pch=24,bg="blue",cex=1.1)}
    }
  }
abline(v=40.5,lty=2, col="darkgrey")

