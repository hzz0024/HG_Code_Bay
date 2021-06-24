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
# Levene is a function to calculate the delta_p with specified selection coefficient. For each s, it will return the delta_p with p ranging from 0-1.
# (-P internal parameter) the vector of allele frequencies in the habitat (i.e. in the larval pool) before selection
# -s the selection coefficient s

Levene <- function(s)
  {
  P <- seq(0,0.99,0.01)
  #P <- seq(0,0.99,0.01)
  #DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s) # original code
  #DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s+1e-10) # add 1e-10 to deal with 0 denominator issue
  DpLev <- s*P*(1-P)/(1 + 2*s*P - s)  
  }


library(ggplot2)

#p <- rep(seq(0,1,0.01),101)
p <- rep(seq(0,0.99,0.01),100)
s <- sort(rep(seq(0,0.99,0.01),100))
#s <- sort(rep(seq(0,0.495,0.005),100)) 
#delta <- unlist(lapply(seq(0,0.495,0.005),Levene))
delta <- unlist(lapply(seq(0,0.99,0.01),Levene))
data <- data.frame(p,s,delta)
write.table(data, "./data.csv", sep=",", quote=F, row.names=T, col.names=NA)
d <- ggplot(data, aes(p, s, z = delta))
d + stat_summary2d() + scale_fill_gradient2(name = "Delta p", high = "red", mid = "grey",low = "black",midpoint = 0.05,guide = "colourbar")


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

#check if odd or even number
is.odd <- function(x) x %% 2 != 0
# function to combine the p-value, multi options are available (this is the new method with correct df estimate)
combinePValues <- function(..., 
                           method=c("fisher", "z", "simes", "berger", "holm-middle"), 
                           weights=NULL, log.p=FALSE, min.prop=0.5)
{
  input <- list(...)
  if (length(input)==1L) {
    return(unname(input[[1]])) # returning directly.
  }
  Np <- unique(lengths(input))
  if (length(Np) != 1) {
    stop("all p-value vectors must have the same length")
  }
  
  method <- match.arg(method)
  switch(method,
         fisher={
           if (log.p) {
             all.logp <- input
           } else {
             all.logp <- lapply(input, FUN=log)
           }
           
           n <- integer(Np)
           X <- numeric(Np)
           for (i in seq_along(all.logp)) {
             current <- all.logp[[i]]
             keep <- !is.na(current)
             X[keep] <- X[keep] + current[keep]
             n <- n + keep
           }
           
           n[n==0] <- NA_real_ # ensure that we get NA outputs.
           pchisq(-2*X, df=2*n, lower.tail=FALSE, log.p=log.p)
         },
         simes=combine_simes(input, log.p),
         `holm-middle`=combine_holm_middle(input, log.p, min.prop),
         z={
           if (is.null(weights)) {
             weights <- rep(1, length(input))
           } else if (length(weights)!=length(input)) {
             stop("'length(weights)' must be equal to number of vectors in '...'")
           } else {
             check <- unlist(lapply(weights, range))
             if (any(is.na(check) | check<=0)) {
               stop("weights must be positive")
             }
           }
           
           Z <- W2 <- numeric(Np)
           for (i in seq_along(input)) {
             current <- input[[i]]
             keep <- !is.na(current)
             Z[keep] <- Z[keep] + qnorm(current[keep], log.p=log.p) * weights[[i]]
             W2 <- W2 + ifelse(keep, weights[[i]]^2, 0)
           }
           
           # Combining p-values of 0 and 1 will yield zscores of -Inf + Inf => NaN.
           # Here, we set them to 0 to get p-values of 0.5, as the Z-method doesn't
           # give coherent answers when you have p-values at contradicting extremes.
           Z[is.nan(Z)] <- 0
           
           W2[W2==0] <- NA_real_ # ensure we get NA outputs.
           pnorm(Z/sqrt(W2), log.p=log.p) 
         },
         berger={
           do.call(pmax, c(input, list(na.rm=TRUE)))
         }
  )
}
############################
# get direction for deltap #
############################

Power <- function(N, K, data)
  {
  d <- data
  d$powerz <- NA
  d$powerfish <- NA
  d$mean_deltap <- NA
  for (r in (1:nrow(data)))
  #for (r in (1:100))
    {
    # first replicate
    RAC1 <- rbinom(100,(2*N),data[r,1])
    RAF1 <- RAC1/(2*N)
    SAC1 <- rbinom(100,(2*K),(rbinom(1,20000,(data[r,1]+data[r,3]))/20000))
    SAF1 <- SAC1/(2*K)
    # second replicate
    RAC2 <- rbinom(100,(2*N),data[r,1])
    RAF2 <- RAC2/(2*N)
    SAC2 <- rbinom(100,(2*K),(rbinom(1,20000,(data[r,1]+data[r,3]))/20000))
    SAF2 <- SAC2/(2*K)
    DAF1 <- SAF1-RAF1
    DAF2 <- SAF2-RAF2
    DP_mean = (DAF1+DAF2)/2
    DP_mean.greater = DP_mean >= 0
    alternative = DP_mean.greater
    alternative[alternative == TRUE] <- 'greater'
    alternative[alternative == FALSE] <- 'less'
    PV <- NULL
    PVfish <- NULL
    PVz <- NULL
    for (t in (1:100))
      {
      if (DP_mean[t]==0)
        {
        PVfish <- c(PVfish, 1)
        }
      if (DP_mean[t]>0)
        {
        #PV <- c(PV, length(which(Ddelta(IAC[t],N,K,theta,10000)>=DAF[t]))/10000)
        CH1 = c(SAC1[t], 2*K-SAC1[t])
        REF1 = c(RAC1[t], (2*N)-RAC1[t])
        M1 = as.table(cbind(CH1, REF1))
        PV1 <- fisher.test(M1, alternative="greater")$p*2 #multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
        CH2 = c(SAC2[t], (2*K)-SAC2[t])
        REF2 = c(RAC2[t], (2*N)-RAC2[t])
        M2 = as.table(cbind(CH2, REF2))
        PV2 <- fisher.test(M2, alternative="greater")$p*2 #multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
        PV1[PV1 > 1] <- 1
        PV2[PV2 > 1] <- 1
        PVfish_cmp = combinePValues(PV1,PV2, method='fisher')
        PVfish <- c(PVfish, PVfish_cmp)
        PVz_cmp = combinePValues(PV1,PV2, method='z') 
        PVz <- c(PVz, PVz_cmp)
        }
      if (DP_mean[t]<0)
        {
        #PV <- c(PV, length(which(Ddelta(IAC[t],N,K,theta,10000)>=DAF[t]))/10000)
        CH1 = c(SAC1[t], 2*K-SAC1[t])
        REF1 = c(RAC1[t], (2*N)-RAC1[t])
        M1 = as.table(cbind(CH1, REF1))
        PV1 <- fisher.test(M1, alternative="less")$p*2 #multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
        CH2 = c(SAC2[t], (2*K)-SAC2[t])
        REF2 = c(RAC2[t], (2*N)-RAC2[t])
        M2 = as.table(cbind(CH2, REF2))
        PV2 <- fisher.test(M2, alternative="less")$p*2 #multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
        PV1[PV1 > 1] <- 1
        PV2[PV2 > 1] <- 1
        PVfish_cmp = combinePValues(PV1,PV2, method='fisher') 
        PVfish <- c(PVfish, PVfish_cmp)
        PVz_cmp = combinePValues(PV1,PV2, method='z') 
        PVz <- c(PVz, PVz_cmp)
        }
      }
    PVfish_adj = p.adjust(PVfish, method = 'BH')
    PVz_adj = p.adjust(PVz, method = 'BH')
    d$powerz[r] <- length(which(PVz_adj<0.05))
    d$powerfish[r] <- length(which(PVfish_adj<0.05))
    d$mean_deltap[r] <- mean(DP_mean)
    print(r)
    }
  d
  }

POW50 <- Power(50,50,data)
write.table(POW50, "./POW50_s0-0.99.csv", sep=",", quote=T, row.names=T, col.names=NA)

# Summary plot (N0=40,N1=100) #
D <- ggplot(POW50, aes(p, s, z = powerz))
D <- D + stat_summary2d() + theme_bw() + labs(x="",y="")
D + scale_fill_gradient2(name="",limits=c(0,100),breaks=c(0,20,40,60,80),high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")

E <- ggplot(POW50, aes(p, s, z = powerfish))
E <- E + stat_summary2d() + theme_bw() + labs(x="",y="")
E + scale_fill_gradient2(name="",limits=c(0,100),breaks=c(0,20,40,60,80),high="red",mid="grey",low="blue",midpoint=40,guide="colourbar")

G <- ggplot(POW50, aes(p, s, z = (powerz-powerfish)))
G <- G + stat_summary2d() + theme_bw() + labs(x="",y="")
G +  scale_fill_gradient2(name="",limits=c(0,20),breaks=c(0,5,10,15,20),high="red",mid="grey",low="blue",midpoint=10,guide="colourbar")


# detailed plots #
D <- ggplot(POW50, aes(p, s, z = powerz))
D + stat_summary2d(bins=100) + theme_bw() + scale_fill_gradient2(name="Power Z",high="red",mid="grey",low="royalblue",midpoint=40,guide="colourbar")
E <- ggplot(POW50, aes(p, s, z = powerfish))
E + stat_summary2d(bins=100) + theme_bw() + scale_fill_gradient2(name="Power Fisher",high="red",mid="grey",low="royalblue",midpoint=40,guide="colourbar")
G <- ggplot(POW50, aes(p, s, z = (powerfish-powerz)))
G + stat_summary2d(bins=100) + theme_bw() + scale_fill_gradient2(name="power difference",high="red",mid="grey",low="royalblue",midpoint=10,guide="colourbar")

###################################### formal plot ###################################### 
install.packages("ggplot2")
library("ggplot2")
library("gridExtra")
require(grid)

file2 = "POW50_s0-0.99.csv"
file2 = "POW_100_50_s0-0.99.csv"
POW50 <- read.delim(file2, header = TRUE, sep=',')

# data delta p
D2 <- ggplot(POW50, aes(p, s, z=mean_deltap))
D2 <- D2 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_0 <- D2+ scale_fill_gradient2(name=expression(Delta~italic(p)),limits=c(-0.1,0.4),breaks=c(-0.1,0.0,0.1,0.2,0.3, 0.4),high="red",mid="grey",low="#B6E5D8",midpoint=0.0,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_0

D2 <- ggplot(POW50, aes(p, s, z = powerz))
D2 <- D2 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_1 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_1

D2 <- ggplot(POW50, aes(s , mean_deltap , z = powerz))
D2 <- D2 + stat_summary2d() + theme_classic() + labs(x="Selection coefficient",y=expression(Delta~italic(p)))
plot50_2 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20))
plot50_2

D2 <- ggplot(POW50, aes(p , mean_deltap , z = powerz))
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



library("ggplot2")
library("gridExtra")
require(grid)
file1 = "POW50_s0-0.99.csv"

POW50 <- read.delim(file1, header = TRUE, sep=',')

D1 <- ggplot(POW50, aes(p, s, z = powerz))
D1 <- D1 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_1 <- D1+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20)) + 
  ggtitle(expression(Combined~Fisher~italic(p)~italic(N)[0]~{textstyle("=")}~50~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50))+
  theme(plot.title = element_text(size=14))
plot50_1

file2 = "POW_100_50_s0-0.99.csv"

POW100_50 <- read.delim(file2, header = TRUE, sep=',')
D2 <- ggplot(POW100_50, aes(p, s, z = powerz))
D2 <- D2 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="")
plot50_2 <- D2+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20)) + 
  ggtitle(expression(Combined~Fisher~italic(p)~italic(N)[0]~{textstyle("=")}~100~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50))+
  theme(plot.title = element_text(size=14))
plot50_2

jpeg("power_comp.jpg", width = 16, height = 9, units = 'in', res = 300)
grid_arrange_shared_legend(plot50_1, plot50_2,ncol=2, nrow = 1, position = "right")
dev.off()

file2 = "POW_200_50_s0-0.99.csv"

POW200_50 <- read.delim(file2, header = TRUE, sep=',')
D3 <- ggplot(POW200_50, aes(p, s, z = powerz))
D3 <- D3 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="Selection coefficient")
plot50_3 <- D3+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20)) + 
  ggtitle(expression(Combined~Fisher~italic(p)~italic(N)[0]~{textstyle("=")}~200~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50))+
  theme(plot.title = element_text(size=14))
plot50_3


file2 = "POW_300_50_s0-0.99.csv"

POW300_50 <- read.delim(file2, header = TRUE, sep=',')
D4 <- ggplot(POW300_50, aes(p, s, z = powerz))
D4 <- D4 + stat_summary2d() + theme_bw() + labs(x=expression(italic(p)*0),y="")
plot50_4 <- D4+ scale_fill_gradient2(name="Power",limits=c(0,1),breaks=c(0,0.2,0.4,0.6, 0.8, 1),high="red",mid="grey",low="#B6E5D8",midpoint=0.5,guide="colourbar") +
  theme(text = element_text(size=20)) + 
  ggtitle(expression(Combined~Fisher~italic(p)~italic(N)[0]~{textstyle("=")}~300~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50))+
  theme(plot.title = element_text(size=14))
plot50_4


jpeg("power_comp.jpg", width = 10, height = 9, units = 'in', res = 300)
grid_arrange_shared_legend(plot50_1, plot50_2, plot50_3, plot50_4,ncol=2, nrow = 2, position = "right")
dev.off()

dat = as.data.frame(cbind(POW50$p, POW50$s, POW50$powerz, POW100_50$powerz))
colnames(dat) <- c("p","s","pow50", "pow100")
G <- ggplot(dat, aes(p, s, z = (pow100-pow50)))
G <- G + stat_summary2d() + theme_bw() + labs(x="",y="")
G + scale_fill_gradient2(name="Power",limits=c(-0.1,0.4),breaks=c(0,0.1,0.2,0.3),high="red",mid="grey",low="#B6E5D8",midpoint=0.1,guide="colourbar")


jpeg("power_diff.jpg", width = 10, height = 9, units = 'in', res = 300)
G + stat_summary2d(bins=100) + theme_bw() + scale_fill_gradient2(name="power",high="red",mid="grey",low="royalblue",midpoint=0.0,guide="colourbar")+
  theme(text = element_text(size=20)) + labs(x=expression(italic(p)*0),y="Selection coefficient")+
  ggtitle(expression(Power~difference~(~italic(N)[0]~{textstyle("=")}~100~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50)~{textstyle("-")}~(italic(N)[0]~{textstyle("=")}~50~{textstyle(",")}~italic(N)[1]~{textstyle("=")}~50)))+
  theme(plot.title = element_text(size=14))
dev.off()

