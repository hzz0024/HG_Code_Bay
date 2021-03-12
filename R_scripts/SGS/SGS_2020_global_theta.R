# load packages
library(gtools)
library(hash)
args = commandArgs(trailingOnly=TRUE)
# to prevent scientific notation caused issue
options(scipen=999)
#######################
# Divide the snp list #
#######################

## the arguments is used to divide the dataset ## 

### PARAMETERS ###

# -1 the ith part of the dataset
# -2 how many parts we want to divide the dataset
# -3 the number of sampling for deltap distribution
# -4 the total number of SNPs

IDX = as.numeric(args[1])
message(paste0('Split:',IDX))

num_split = as.numeric(args[2])
if(is.null(args[2])){
  num_split = 100
}
message(paste0('num_split:',num_split))

nreps = as.numeric(args[3])
if(is.null(args[3])){
  nreps = 10000
}
message(paste0('nreps:',nreps))

total_n = as.numeric(args[4])
if(is.null(args[4])){
  total_n = 2274254
}
message(paste0('total number of locus:',total_n))

interval = total_n/num_split
idxs = seq(interval*IDX+1, interval*(IDX+1))
message(idxs[1])
message(idxs[length(idxs)])

# temp settings
total_n = 2274254

#######################
#   LOAD MAF FILES    #
#######################
# load reference file with header in it
ref = 'REF19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.mafs'
dat_ref <- read.delim(ref, header = TRUE, sep='\t')
ref_n = dat_ref$nInd
ref_k = round(dat_ref$knownEM*dat_ref$nInd*2)

# load challenge file with header in it
chr = 'CHR19_maf0.05_minq20_minmq30_pctind0.7_CV30_masked.mafs'
dat_ch <- read.delim(chr, header = TRUE, sep='\t')
ch_n <- dat_ch$nInd
ch_k <- round(dat_ch$knownEM*dat_ch$nInd*2)

# load window sized theta file with no header in it
#pi_file = 'pi_correct_all.txt'
#pi_dat <- read.delim(pi_file, header = FALSE, sep='\t') # from ch - ref

outputfile = paste0('p_values_theta', IDX, '.txt')
message(outputfile)
#message(dim(dat_ref)[1])
message("Total number of SNP is ", total_n)
#dat_ref = dat_ref[dat_ref$chromo==5,]
#dat_ch= dat_ch[dat_ch$chromo==5,]

######################
# AFS AT EQUILIBRIUM #
######################

## Now use equation 50 in Tajima 1989 to get the Allele        ##
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

Dfreq <- function(X,N,theta){
  theta = as.numeric(theta)
  # p is the vector that stores all the frequencies that the
  # allele can take in a sample of size N
  p=(1:(2*N-1))/(2*N)
  # we get P(X|p) over all possible values of p within the pool
  # this is equivalent to exp(lfactorial(N)-(lfactorial(N - X)+lfactorial(X))) * p^X * (1-p)^(N-K)
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

##########################
#       Theta map        #
##########################

## Save the window sized theta   ##

#theta_map = hash()

#save_theta <- function(pi_dat){
#  for(row in 1:nrow(pi_dat)){
#    
#    if(is.na(pi_dat[row, "V6"])){
#      next
#    }
#    cur_chorm = pi_dat[row, "V1"]
#    cur_pos = pi_dat[row, "V2"]
#    cur_theta = formatC(pi_dat[row, "V6"], digits = 5, format = "f")
#    theta_map[[paste0(cur_chorm,' ',cur_pos)]] = cur_theta
#  }
#}

#save_theta(pi_dat)

################################
# PROBABILITY DENSITY OF DELTA #
################################

## Now randomly draw allele frequency values from    ##
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
  delta = samplefreq - (X/(2*N)) 
}

################################
#    QUANTILE TEST FUNCTION    #
################################

##                function to perform the quantile test            ##
##https://people.stat.sc.edu/hitchcock/Rexamples518section3_2.txt  ##
##            https://freakonometrics.hypotheses.org/4199          ##


quantile.test<-function(x,xstar=0,quantile=.5,alternative="two.sided"){
  n<-length(x)
  p<-quantile
  T1<-sum(x<=xstar)
  T2<-sum(x< xstar)
  if (alternative=="quantile.less") {
    p.value<-1-pbinom(T2-1,n,p)}
  if (alternative=="quantile.greater"){
    p.value<-pbinom(T1,n,p)}
  if (alternative=="two.sided"){
    p.value<-2*min(1-pbinom(T2-1,n,p),pbinom(T1,n,p))}
  list(xstar=xstar,alternative=alternative,T1=T1,T2=T2,p.value=p.value)}

################################
#    OUTLIER IDENTIFICATION    #
################################

# main content of the function

del_cnt = 0
dic <- hash()
p_values = c()
sink(outputfile)
idxs = seq(interval*IDX+1, interval*(IDX+1))
#idxs = seq(1,dim(dat_ref)[1])
idxs = seq(1,500) #test purpose
cnt = 0
for(i in idxs){
  # print the running process
  s = paste0(i,'/',dim(dat_ref)[1])
  message(s,"\r",appendLF=FALSE)
  
  ref_n = dat_ref$nInd[i]
  ref_k = round(dat_ref$knownEM[i]*dat_ref$nInd[i]*2)
  ch_n = dat_ch$nInd[i]
  
  #theta = theta_map[[theta_key]]
  theta = 0.015965112 # replace the local theta with global theta (i.e. 2019 reference theta)
  
  obs_delta=dat_ch$knownEM[i] - dat_ref$knownEM[i]
  delta_ps = Ddelta(X=ref_k, N=ref_n, K=ch_n, theta=theta, nreps=nreps)
  
  # quantile test for p-value
  if(obs_delta>0){
    res = quantile.test(delta_ps,xstar=obs_delta, quantile=0.95, alternative="quantile.less")
    p_value=res$p.value
  }else if(obs_delta<0){
    res = quantile.test(delta_ps,xstar=obs_delta, quantile=0.05, alternative="quantile.greater")
    p_value=res$p.value
  }else{
    message('obs_delta is 0')
    p_value=1.0
  }

  cat(dat_ref$chromo[i])
  cat('\t')
  cat(dat_ref$position[i])
  cat('\t')
  cat(p_value)
  cat('\t')
  
  if(p_value<0.05){
    hist(delta_ps, breaks=30, xlab=paste0("delta_p: " ,dat_ref$chromo[i], "_", dat_ref$position[i]))
    abline(v=obs_delta)
    
    #message(obs_delta)
    cnt = cnt + 1
    cat(1)
    
  }else{
    cat(0)
  }
  cat('\n')
  flush.console()
  
  p_values = c(p_values, p_value)
}
sink()
message(cnt)

#hist(p_values)

#library("qvalue")
#qobj <- qvalue(p = p_values,pi0 = 1)
#lfdr <- qobj$lfdr
#length(lfdr[lfdr<0.01])
#p_values <- read.delim(outputfile,header=FALSE) #read p values from p_values.txt
#out = data.frame(chromo=dat_ref$chromo, position=dat_ref$position, p_value=p_values)
#write.table(out, file = "p_value_list_all.txt", sep = "\t", row.names = FALSE, col.names = FALSE)



