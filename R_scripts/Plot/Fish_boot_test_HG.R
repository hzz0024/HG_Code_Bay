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

do.plot = function(pval.df, MAIN, LEFT, Y, LEN, threshold=0.1){##NOTE THIS FUNCTION IS REPEATED BELOW
  p = get.combined(pval.df)
  DF = data.frame(VARIANT, LG, CM, p)
  
  CEX = 1.2
  adj = p.adjust(DF$p, method = 'BH')
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
  plot(-log(adj, 10)~CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Chromosome", ylab = expression('-log'['10']*'(fdr)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  par(lend = 2)
  if (PLOT.SCALE == TRUE){
    segments(LEFT, Y, (LEFT + LEN), Y, lwd = 2)
    TEXT = paste(LEN, "cM")
    labPos = LEFT + .5 * LEN
    #text(labPos, y = (Y + .4), labels = TEXT)
  }
  rownames(DF) = DF$VARIANT
  #return(DF)
}

## Produce datasets for manhattan plots 

dat = read.delim(file1, header=FALSE, sep='\t')
ps <- combine_p(file1, file2)

# decide whether use test sample or not
if(SAMPLE==TRUE){
  #idx = seq(1,50000) # sampling first 50000
  idx = sample(seq(1,length(dat[,1])), 500000)
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

reverse_pos <- function(genpos){
  chr = sum(genpos > posshift)
  pos = genpos - posshift[chr]
  return(paste0(chr,'_',pos))
}

LG = as.numeric(factor(dat$V1))
CM = genpos
VARIANT = paste0(dat$V1,'_',dat$V2)
lg.names = 1:10
point.size = .3
setYlim = c(-.2, 8)
p = get.combined(ps)
DF = data.frame(VARIANT, LG, CM, p)

#SET UP LABELS FOR THE 10 LINKAGE GROUPS
lgs = 1:10

#GATHER THE ENDPOINTS AND MIDPOINTS OF EACH LINKAGE GROUP FOR PLOTTING PURPOSES
endpoints = c(0)
for (i in lgs){
  sub = CM[LG == i]
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
##################################################################################
################# DO BOOTSTRAP TO DEMONSTRATE HOW RARE PEAKS ARE #################
##################################################################################

##step1. get the null distributions for counts of markers below a p theshold
##note that the p  values in the dataframe are not transformed (i.e. not fdr corrected), and these are used to get the window p values
reps = 10^5
grab.size = 15
cut = 0.05

get.bootstrap = function(DF, id, grab.size, cut){
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
#Divides the linkage groups into windows each with  number of markers equal to grab.size
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
plot.boot.bars = function(DF, DF.ps, grab.size, cut, null.dist, MAIN, bar.Ycoord, bar.width, easy.cut, hard.cut){
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
  for (i in 1:10){
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
  result1$adj.ps = p.adjust(result1$boot.ps, method = "BH")
  easy = result1[result1$adj.ps > hard.cut,]
  easy = easy[easy$adj.ps < easy.cut,]
  hard = result1[result1$adj.ps < hard.cut,]
  print(head(easy))
  print(nrow(easy))
  do.plot(DF.ps, MAIN, 100, 5, 100, threshold = 0.05)
  par(xpd=TRUE)
  top = bar.Ycoord + 0.5*bar.width
  bottom = bar.Ycoord - 0.5*bar.width
  print(seq(bottom, top, by = 0.01))
  for (y in seq(bottom, top, by = 0.01)){
    if (nrow(easy) >= 1 && FALSE){
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
    
    return(output)
  }
}


ps.boot.df = get.bootstrap(DF, 'ps_sample', grab.size, cut)
hist(ps.boot.df$ps)

ps.boot = ps.boot.df$counts

#######################################################
#################### DO FINAL PLOT ####################
#######################################################
###SET UP SOME PLOTTING VARIABLES COLOR CODING SCHEME
quartz()
par(mfrow = c(1,1))
EASY.BAR.COLOR = 'darkgreen'
EASY.CUT = 0.05
HARD.BAR.COLOR = 'green'
HARD.CUT = 0.01
PLOT.SCALE = TRUE
#PLOT
setYlim = c(-.75, 4)
point.size = .75
par(mar=c(2,2,1,2)+2)

#png(paste0(GLOBAL_NAME, 'bootstrap.jpg'))
#bars = plot.boot.bars(DF, ps, grab.size, cut, ps.boot, 'Test', .825*setYlim[1], -.75*setYlim[1], EASY.CUT, HARD.CUT)
#dev.off()

#lefts = lapply(bars$lefts, reverse_pos)
#rights = lapply(bars$rights, reverse_pos)
#outputs = cbind(lefts, rights)
#write.table(outputs, paste0(GLOBAL_NAME,'output_bars.txt'))

png(paste0(GLOBAL_NAME, 'bootstrap1.jpg'))
bars1 = plot.boot.bars(DF, ps, grab.size, cut, ps.boot, 'Test', .825*setYlim[1], -.75*setYlim[1], 0.1, 0.05)
dev.off()

lefts = sapply(bars1$lefts, reverse_pos)
rights = sapply(bars1$rights, reverse_pos)
outputs = data.frame(lefts, rights, bars1$cut)
colnames(outputs) = c('left', 'right', 'cut')
write.table(outputs, paste0(GLOBAL_NAME,'output_bars1.txt'))

#png(paste0(GLOBAL_NAME, 'bootstrap2.jpg'))
#bars2 = plot.boot.bars(DF, ps, grab.size, cut, ps.boot, 'Test', .825*setYlim[1], -.75*setYlim[1], 0.05, 0.01)
#dev.off()

#lefts = sapply(bars2$lefts, reverse_pos)
#rights = sapply(bars2$rights, reverse_pos)
#outputs = data.frame(lefts, rights, bars2$cut)
#colnames(outputs) = c('left', 'right', 'cut')
#write.table(outputs, paste0(GLOBAL_NAME,'output_bars2.txt'))


