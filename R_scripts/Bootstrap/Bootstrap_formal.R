##########################################################################
########################### FISHER'S TESTS ###############################
##########################################################################
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
    combined.p = combinePValues(x,y, method='fisher')
    #combined.p = combinePValues(x,y, method='z')
    combined.p.values = append(combined.p.values, combined.p)
  }
  return(combined.p.values)
}

# Combined p-value adjustment and output the SNPs with fdr < 0.01 and < 0.05 @HG
# ARGS: ps - a vector of combined p-values; dat - file with exact test p-values; out1/out5 - output names for SNPs with fdr < 0.01 and fdr < 0.05
#fisher_method <- function(ps, dat, outall, out1, out5){
#  
#  ps = get.combined(ps)
#  adj = p.adjust(ps, method = 'BH')
#  
# #idx = adj
#  datall = paste0(dat$V1,'\t',dat$V2,'\t',adj)
#  message(paste0('all: ',length(dat$V1)))
#  write.table(datall, outall, quote = FALSE, row.names = FALSE)
#  
#  idx = adj < 0.05
#  dat5 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
#  message(paste0('0.05: ',length(dat$V1[idx])))
#  write.table(dat5, out5, quote = FALSE, row.names = FALSE)
#  
#  idx = adj < 0.01
#  dat1 = paste0(dat$V1[idx],'\t',dat$V2[idx],'\t',adj[idx])
#  message(paste0('0.01: ',length(dat$V1[idx])))
#  write.table(dat1, out1, quote = FALSE, row.names = FALSE)
#}

do.plot = function(DF, MAIN, threshold=0.05){##NOTE THIS FUNCTION IS REPEATED BELOW
  adj = DF$p #adjusted conbined p
  CEX = 1.2
  
  #set up grey and black points alternating by LG
  colors = is.odd(DF$LG)
  colors[colors == TRUE] <- 'black'
  colors[colors == FALSE] <- 'grey'
  #make significant points red
  x = adj < threshold
  for (i in 1:length(x)){
    if (x[i] == TRUE){
      colors[i] <- 'red'
    }
  }
  # comment out if all points need to be plotted
  #DF = DF[x==TRUE,] #@HG
  #colors = colors[x==TRUE] #@HG
  #CM = CM[x==TRUE] #@HG
  
  
  
  #GATHER THE ENDPOINTS AND MIDPOINTS OF EACH LINKAGE GROUP FOR PLOTTING PURPOSES
  endpoints = c(0)
  for (i in lgs){
    sub = DF$CM[DF$LG == i]
    end = max(sub)
    endpoints = append(endpoints, end)
  }
  mids = c()
  for (i in seq(1:(length(endpoints) - 1))){
    left = endpoints[i]
    right = endpoints[(i + 1)]
    point = (left + right) / 2
    mids = append(mids, point)
  }
  
  plot(-log(adj, 10)~DF$CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(fdr)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  rownames(DF) = DF$VARIANT
}



##################################################################################
################# DO BOOTSTRAP TO DEMONSTRATE HOW RARE PEAKS ARE #################
##################################################################################
get.bootstrap = function(DF, id, grab.size, cut, reps = 10^5){
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
  write.table(results, paste(GLOBAL_NAME, id, grab.size, cut, sep = "_"), row.names = F, quote = F)
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

#FUNCTION setup.windows
#Divides the linkage groups into windows each with number of markers equal to grab.size
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

#FUNCTION plot.boot.bars
#This combines the do.plot() function used above with the bootstrapping functions to plot
#a finalized figure that shows the combined p values and color coded bars below regions
#with significant bootstrapping signal. Also outputs the bootstrap coordiantes for later use.
plot.boot.bars = function(DF, grab.size, cut, null.dist, MAIN, bar.Ycoord, bar.width, easy.cut, hard.cut){
  #ARGUMENTS: 
  # DF = the dataframe with the p values output for a particular cross (output from function do.plot()) 
  # DF.ps = data frame with the one-sided p values for each individual replicate
  # grab.size = the window size (in number of loci) used to build the null distribution beging used
  # cut = the cutoff p value cutoff for the null distribution
  # null.dist = the null distribution output from get.bootstrap()
  # MAIN = title
  # barYcoord = the vertical coordinate for placing the bootstrap bars
  # bar.width = the vertical width of the bootstrap bar
  # easy.cut = the unconservative p value cutoff for bootstrap bars
  # hard.cut = the conservative cutoff
  x = c()
  boot.ps = c()
  sig.counts = c()
  window.left.bounds = c()
  window.right.bounds = c()
  for (i in lgs){
    lg.sub = DF[DF$LG == i,] ##pull out the linkage groups one at a time
    ##set up the windows. This is more complicated than just seq(by = grab.size, because we want to include every marker)
    window.sides = setup.windows(grab.size, length(lg.sub[,1]))##sets up the windows as a dataframe of left and right sides
    for (win in seq(1:(length(window.sides[,1])))){
      left = window.sides[win,1]
      right = window.sides[win,2]
      sub = lg.sub[left:right,]
      win.left.cM = sub[1,'CM']##grabs the left side of this window in cM
      win.right.cM = sub[length(sub[,1]),'CM']##grabs the right side of the window
      window.left.bounds = append(window.left.bounds, win.left.cM)
      window.right.bounds = append(window.right.bounds, win.right.cM)
      #get p value for counts of significant markers (turn this on and off)
      sig.sub = sub$p[sub$p < cut]
      obs.count = length(sig.sub)
      p = get.boot.p(null.dist, obs.count)
      boot.ps = append(boot.ps, p)
      cM = median(sub$CM)
      x = append(x, cM)
    }
  }
  result1 = data.frame(x, boot.ps, window.left.bounds, window.right.bounds)
  #result1$adj.ps = p.adjust(result1$boot.ps, method = "BH")
  result1$adj.ps = result1$boot.ps
  #write.table(result1, paste(GLOBAL_NAME, grab.size, cut, sep = "_"), row.names = F, quote = F)
  easy = result1[result1$adj.ps > hard.cut,]
  easy = easy[easy$adj.ps < easy.cut,]
  hard = result1[result1$adj.ps < hard.cut,]
  print(head(easy))
  print(nrow(easy))
  do.plot(DF, MAIN, threshold = 0.05)
  
  
  par(xpd=TRUE)
  top = bar.Ycoord + 0.5*bar.width
  bottom = bar.Ycoord - 0.5*bar.width
  print(seq(bottom, top, by = 0.01))
  for (y in seq(bottom, top, by = 0.01)){
    if (nrow(easy) >= 1){
      for (i in 1:nrow(easy)){
        segments(easy$window.left.bounds[i], y, easy$window.right.bounds[i], y, col = EASY.BAR.COLOR, lwd = 2)
      }
    }
    if (nrow(hard) >= 1){
      for (i in 1:nrow(hard)){
        segments(hard$window.left.bounds[i], y, hard$window.right.bounds[i], y, col = HARD.BAR.COLOR, lwd = 2)
      }
    }
  }
  par(xpd=FALSE)
  if (nrow(hard) > 0 || nrow(easy) > 0){
      print("easy:")
      print(easy)
      print("hard:")
      print(hard)
      COLNAMES = c('lefts', 'rights', 'cut')
      
      if(nrow(easy) > 0){
        e.lefts = easy$window.left.bounds
        e.rights = easy$window.right.bounds
        output.e = data.frame(e.lefts, e.rights)
        output.e$cut = "easy"
        if(nrow(hard) == 0){
          colnames(output.e) = COLNAMES
          output = output.e
        }
      }
      
      if(nrow(hard) > 0){
        h.lefts = hard$window.left.bounds
        h.rights = hard$window.right.bounds
        output.h = data.frame(h.lefts, h.rights)
        output.h$cut = "hard"
        if(nrow(easy) == 0){
          colnames(output.h) = COLNAMES
          output = output.h
        }
      }
      
      if(nrow(hard) > 0 && nrow(easy) > 0){
        colnames(output.e) = COLNAMES
        colnames(output.h) = COLNAMES
        output = rbind(output.e, output.h)
      }
      
      #return(output)
  }
  #return(result1)
  return(list(result1, output))
}

##########################################################################
########################### MAIN FUNCTION ################################
##########################################################################
##############--PARAMETER--##############
SAMPLE = FALSE
GLOBAL_NAME = 'REF-CH-SR-HC_'
file1 = "fish_REF_CH_REF-CH-SR-HC.txt"
file2 = "fish_SR_HC_REF-CH-SR-HC.txt"
#############--MAIN--###################
dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

if(SAMPLE==TRUE){
  idx = seq(1200000,1300000) 
  ps = ps[idx,]
  rownames(ps) <- seq(length=nrow(ps))
  dat = dat[idx,]
  rownames(dat) <- seq(length=nrow(dat))
}

dat$V1 = as.numeric(factor(dat$V1))
idxs = order(dat$V1, dat$V2)
dat = dat[idxs,]
ps = ps[idxs,]

pos = dat$V2
chr = dat$V1
posmin <- tapply(pos,chr, min);
posmax <- tapply(pos,chr, max);
POSSHIFT <- head(c(0,cumsum(posmax)),-1);
names(POSSHIFT) <- levels(chr)
genpos <- pos + POSSHIFT[chr];

reverse_pos <- function(genpos){
  chr = sum(genpos > POSSHIFT)
  pos = genpos - POSSHIFT[chr]
  return(paste0(chr,'_',pos))
}


#LG = as.numeric(factor(dat$V1))
LG = dat$V1
CM = genpos
VARIANT = paste0(dat$V1,'_',dat$V2)
lgs = 1:10
lgs = min(LG):max(LG)
lg.names = lgs
print(lgs)
point.size = .3
setYlim = c(-.2, 8)
p = get.combined(ps)
#@change: adjust p
p = p.adjust(p, 'BH')
DF = data.frame(VARIANT, LG, CM, p)


cut = 0.05
EASY.BAR.COLOR = 'green'
EASY.CUT = 0.1
HARD.BAR.COLOR = 'darkgreen'
HARD.CUT = 0.05
for(grab.size in c(15, 25, 50)){
  
  ps.boot.df = get.bootstrap(DF, 'ps_sample', grab.size, cut, reps = 10^5)
  
  png(paste0(GLOBAL_NAME, grab.size, '_bootstrap.jpg'))
  results = plot.boot.bars(DF, grab.size, cut, ps.boot.df$counts, paste0(GLOBAL_NAME,'barplot') , .825*setYlim[1], -.75*setYlim[1], EASY.CUT, HARD.CUT)
  dev.off()
  
  result = results[1]
  result = data.frame(result)
  write.csv(result, paste0(GLOBAL_NAME, grab.size, '_grab.size.csv'))
  
  bars1 = results[2]
  bars1 = data.frame(bars1)
  lefts = sapply(bars1$lefts, reverse_pos)
  rights = sapply(bars1$rights, reverse_pos)
  outputs = data.frame(lefts, rights, bars1$cut)
  colnames(outputs) = c('left', 'right', 'cut')
  write.table(outputs, paste0(GLOBAL_NAME, grab.size,'_autput_bars.txt'))
}

