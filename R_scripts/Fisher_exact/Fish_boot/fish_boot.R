SAMPLE = TRUE
GLOBAL_NAME = 'hard_SR-REF-COH-ARN_'
file1 = "fish_SR_REF_SR-REF-COH-ARN.txt"
file2 = "fish_COH_ARN_SR-REF-COH-ARN.txt"
#================================

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
    # combined.p = combinePValues(x,y, method='fisher')
    combined.p = combinePValues(x,y, method='z')
    combined.p.values = append(combined.p.values, combined.p)
  }
  return(combined.p.values)
}

# Functio to combine two p-values @HG
# ARGS: p1_name/p2_name - p-values from fish exact test
combine_p <- function(p1_name, p2_name){
  dat_1 <- read.delim(p1_name, header = FALSE, sep='\t')
  dat_2 <- read.delim(p2_name, header = FALSE, sep='\t')
  p1 <- dat_1$V7
  p2 <- dat_2$V7
  ps <- data.frame(p1, p2)
  return(ps)
}

setup.windows = function(grab.size, lg.tot){
  first.whole = floor(lg.tot/grab.size)##sets up the first windows up until the last
  left.sides = c()
  right.sides = c()
  left = 0
  for (i in 1:first.whole){
    right = grab.size*i
    left.sides = append(left.sides, left)
    right.sides = append(right.sides, right)
    left = right  ##now make the new left the previous right side, this way all adjacent windows will overlap
  }
  last.right = lg.tot
  last.left = lg.tot - (grab.size-1)
  left.sides = append(left.sides, last.left)
  right.sides = append(right.sides, last.right)
  window.sides = data.frame(left.sides, right.sides)
  return(window.sides)
}

## Produce datasets for manhattan plots 
dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

# decide whether use test sample or not
if(SAMPLE==TRUE){
  #idx = seq(1,50000) # sampling first 50000
  idx = sample(seq(1,length(dat[,1])), 50000)
  ps = ps[idx,]
  rownames(ps) <- seq(length=nrow(ps))
  dat = dat[idx,]
  rownames(dat) <- seq(length=nrow(dat))
}
idxs = order(dat$V1, dat$V2)
dat = dat[idxs,]
ps = ps[idxs,]

pos = dat$V2
chr = as.numeric(factor(dat$V1))
posmin <- tapply(pos,chr, min);
posmax <- tapply(pos,chr, max);
posshift <- head(c(0,cumsum(posmax)),-1);
names(posshift) <- levels(chr)
genpos <- pos + posshift[chr];

LG = as.numeric(factor(dat$V1))
CM = genpos
VARIANT = paste0(dat$V1,'_',dat$V2)
lg.names = 1:10
point.size = .3
setYlim = c(-.2, 8)

p = get.combined(ps)
DF = data.frame(VARIANT, LG, CM, p)


##################################################################################
################# DO BOOTSTRAP TO DEMONSTRATE HOW RARE PEAKS ARE #################
##################################################################################

##step1. get the null distributions for counts of markers below a p theshold
##note that the p  values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values
reps = 10^5
grab.size = 15
cut = 0.05

get.bootstrap = function(DF, grab.size, cut){
  counts = c()
  ps = c()
  for (i in 1:reps){
    sub = sample(DF$p, grab.size)
    sub2 = sub[sub < cut]
    new.count = length(sub2)
    new.p = mean(sub)
    counts = append(counts, new.count)
    ps = append(ps, new.p)
  }
  results = data.frame(counts, ps)
  #write.table(results, paste(GLOBAL_NAME, id, grab.size, cut, sep = "_"), row.names = F, quote = F)
  return(results)
}##Arugments: DF = the dataframe; id = a string to id the file from when exported

#SET UP THE FUNCTIONS FOR PERFORMING THE BOOTSTRAPPING

# FUNCTION get.boot.p
# A putatively selected variant is one with an unadjustd p value 0.05
# This function compares a count of putatively selected variants
# To the null distribution and returns the proportion of null entires 
# with equal more greater number of putatively selected variants
get.boot.p = function(null, obs){
  tot = length(null)
  more.extreme.count = length(null[null >= obs])
  p.value = more.extreme.count/tot
  return(p.value)
}

ps.boot.df = get.bootstrap(DF, grab.size, cut)
hist(ps.boot.df$ps)

null.dist = ps.boot.df$counts
x = c()
boot.ps = c()
sig.counts = c()
window.left.bounds = c()
window.right.bounds = c()
for (i in 1:10){
  lg.sub = DF[DF$LG == i,] ##pull out the linkage groups one at a time
  ##set up the windows. This is more complicated than just seq(by = grab.size, because we want to include every marker)
  window.sides = setup.windows(grab.size, length(lg.sub[,1]))##sets up the windows as a dataframe of left and right sides
  for (win in seq(1:(length(window.sides[,1])))){
    left = window.sides[win,1]
    right = window.sides[win,2]
    print(left)
    print(right)
    sub = lg.sub[left:right,]
    #get p value for counts of significant markers (turn this on and off)
    sig.sub = sub$p[sub$p < cut]
    obs.count = length(sig.sub)
    p = get.boot.p(null.dist, obs.count)
    boot.ps = append(boot.ps, p)
  }
}

adj.ps = p.adjust(boot.ps, method = "BH")
sum(adj.ps < 0.05)


