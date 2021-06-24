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
library(ggplot2)
Levene <- function(s)
{
  P <- seq(0,0.99,0.01)
  #P <- seq(0,0.99,0.01)
  #DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s) # original code
  #DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s+1e-10) # add 1e-10 to deal with 0 denominator issue
  DpLev <- s*P*(1-P)/(1 + 2*s*P - s)  
}

p <- seq(0,0.99,0.01)
DPLev <- Levene(0.1)
df <- data.frame(cbind(p,DPLev))

d <- ggplot(df, aes(x=p, y=DPLev))
f1 <- d + geom_point() + theme_bw() + labs(x=expression(italic(p)*0),y=expression(Levene~Delta~italic(p))) +
  ggtitle("Selection coefficient = 0.1")

DPLev <- Levene(0.2)
df <- data.frame(cbind(p,DPLev))

d <- ggplot(df, aes(x=p, y=DPLev))
f2 <- d + geom_point() + theme_bw() + labs(x=expression(italic(p)*0),y=expression(Levene~Delta~italic(p))) +
  ggtitle("Selection coefficient = 0.2")

DPLev <- Levene(0.5)
df <- data.frame(cbind(p,DPLev))

d <- ggplot(df, aes(x=p, y=DPLev))
f3 <-  d + geom_point() + theme_bw() + labs(x=expression(italic(p)*0),y=expression(Levene~Delta~italic(p))) +
  ggtitle("Selection coefficient = 0.5")

DPLev <- Levene(0.9)
df <- data.frame(cbind(p,DPLev))

d <- ggplot(df, aes(x=p, y=DPLev))
f4 <- d + geom_point() + theme_bw() + labs(x=expression(italic(p)*0),y=expression(Levene~Delta~italic(p))) +
  ggtitle("Selection coefficient = 0.9")

library("ggplot2")
library("gridExtra")
require(grid)
grid.arrange(f1, f2, f3, f4, nrow = 2)


