SAMPLE = FALSE
GLOBAL_NAME = 'REF-CH-SR-HC_'
file1 = "fish_REF_CH_REF-CH-SR-HC.txt"
file2 = "fish_SR_HC_REF-CH-SR-HC.txt"
#================================

#performs fisher's method to combine a set of p values (this is old method with an error in the df part)
#ARGUMENTS: ps = a vector of p values you want combined 
fisher.method=function(ps) {
  ftest=-2*sum(log(ps))
  #df=length(2*ps)
  df=2*length(ps)
  pv=1-pchisq(q=ftest,df=df)
}

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

#FUNCTION get.combined():
#runs fisher.method() on a dataframe of p values
#ARGUMENTS: DF = a two column data frame with p values in it
get.combined = function(DF){
  combined.p.values = c()
  alt = c()
  for (i in seq(1, length(DF[,1]))){ #(apparently did not know of nrow() yet)
    x = DF[i,1]
    y = DF[i,2]
    #combined.p = fisher.method(c(x,y))
    #combined.p = combinePValues(x,y, method='fisher')
    combined.p = combinePValues(x,y, method='z')
    combined.p.values = append(combined.p.values, combined.p)
  }
  return(combined.p.values)
}

# Functio to combine two p-values @HG
# ARGS: p1_name/p2_name - files with exact test p-values
combine_p <- function(p1_name, p2_name){
  dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
  dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
  p1 <- dat_1$V7
  p2 <- dat_2$V7
  ps <- data.frame(p1, p2)
  return(ps)
}

# Combined p-value adjustment and output the SNPs with fdr < 0.01 and < 0.05 @HG
# ARGS: ps - a vector of combined p-values; dat - file with exact test p-values; out1/out5 - output names for SNPs with fdr < 0.01 and fdr < 0.05
fisher_method <- function(ps, dat, outall, out1, out5){
  
  ps = get.combined(ps)
  adj = p.adjust(ps, method = 'BH')
  
  #idx = adj
  datall = paste0(dat$V1,'\t',dat$V2,'\t',adj)
  message(paste0('all: ',length(dat$V1)))
  write.table(datall, outall, quote = FALSE, row.names = FALSE)
  
  idx = adj < 0.05
  dat5 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
  message(paste0('0.05: ',length(dat$V1[idx])))
  write.table(dat5, out5, quote = FALSE, row.names = FALSE)
  
  idx = adj < 0.01
  dat1 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
  message(paste0('0.01: ',length(dat$V1[idx])))
  write.table(dat1, out1, quote = FALSE, row.names = FALSE)
}

## Produce datasets for manhattan plots 

dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

# decide whether use test sample or not
if(SAMPLE==TRUE){
  idx = seq(1,10000) #del!!!
  #idx = sample(seq(1,length(dat[,1])), 50000)
  ps = ps[idx,]
  rownames(ps) <- seq(length=nrow(ps))
  dat = dat[idx,]
  rownames(dat) <- seq(length=nrow(dat))
}

idxs = order(dat$V1, dat$V2)
dat = dat[idxs,]
ps = ps[idxs,]

fisher_method(ps, dat, outall=paste0(GLOBAL_NAME, 'out_all_z.txt'), out1=paste0(GLOBAL_NAME, 'out_0.01_z.txt'), out5=paste0(GLOBAL_NAME, 'out_0.05_z.txt')) #write outliers to output files, no return
